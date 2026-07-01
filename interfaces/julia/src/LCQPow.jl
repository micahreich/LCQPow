module LCQPow

using CxxWrap
using Preferences
using LinearAlgebra: dot

# Path to the compiled CxxWrap library. Override with a LocalPreferences.toml
# entry `liblcqpow_julia_path` if you build somewhere other than `build/lib`.
const _liblcqpow_julia = @load_preference("liblcqpow_julia_path",
    joinpath(@__DIR__, "..", "..", "..", "build", "lib", "liblcqpow_julia"))

@wrapmodule(() -> _liblcqpow_julia)

function __init__()
    @initcxx
end

# Public enum-style constants (mirror LCQPow::QPSolver / LCQPow::PrintLevel).
const QPOASES_DENSE  = 0
const QPOASES_SPARSE = 1
const OSQP_SPARSE    = 2

const PRINT_NONE  = 0
const PRINT_OUTER = 1
const PRINT_INNER = 2

# LCQPow treats bounds with magnitude >= this as unbounded (LCQPow::Utilities::INFTY).
const INFTY = 1.0e20

# Solver options, in the exact order expected by the C++ binding. Values are
# passed as a Float64 vector where NaN means "use the LCQPow default".
const OPTION_NAMES = (
    :stationarity_tolerance,
    :complementarity_tolerance,
    :initial_penalty_parameter,
    :penalty_update_factor,
    :solve_zero_penalty_first,
    :perturb_step,
    :max_iterations,
    :max_penalty_parameter,
    :n_dynamic_penalty,
    :eta_dynamic_penalty,
    :print_level,
    :store_steps,
    :qp_solver,
)

const OPTION_INDEX = Dict(name => i for (i, name) in pairs(OPTION_NAMES))

# LCQPow's status/message strings are formatted for console printing (stray '#'
# and trailing newlines); tidy them up for programmatic use.
_clean(s) = strip(replace(String(s), '#' => "", '\n' => ' '))

function _option_values(kwargs)
    values = fill(NaN, length(OPTION_NAMES))
    for (name, value) in pairs(kwargs)
        index = get(OPTION_INDEX, name, nothing)
        isnothing(index) && throw(ArgumentError("unknown LCQPow option: $name"))
        values[index] = Float64(value)
    end
    return values
end

"""
    solve_qpcc_with_lcqpow(data::NamedTuple; x0=nothing, print_level=PRINT_NONE, use_dummy_constant=false, kwargs...) -> NamedTuple

Solve Marble-form QPCC data (the CRISP formulation) with LCQPow:

    minimize    1/2 x'Q x + q'x + c0
    subject to  J_eq x + b_eq == 0
                J_ineq x + b_ineq >= 0
                0 <= L x + l   complements   R x + r >= 0

`data` is a NamedTuple with fields `Q`, `q`, `c0`, `J_eq`, `b_eq`, `J_ineq`,
`b_ineq`, `L`, `l`, `R`, `r`. Pass zero-row blocks (e.g. `zeros(0, n)`) for
constraint types you don't have. This maps onto LCQPow's native form as:

  - equalities   `J_eq x == -b_eq`     -> two-sided linear row with `lbA = ubA = -b_eq`;
  - inequalities `J_ineq x >= -b_ineq` -> linear row with `lbA = -b_ineq, ubA = +∞`;
  - complementarity -> `lbL = -l`, `lbR = -r` (LCQPow enforces `(Lx-lbL) ⊥ (Rx-lbR)`).

`x0` is an optional primal warm start. Choose the QP subsolver with
`qp_solver = LCQPow.QPOASES_DENSE | QPOASES_SPARSE | OSQP_SPARSE`. Other solver
hyperparameters are forwarded as keyword arguments (see `OPTION_NAMES`), e.g.
`stationarity_tolerance`, `max_iterations`.

LCQPow can be unreliable when the complementarity constants `l`, `r` are nonzero
(they become nonzero `lbL`/`lbR` lower bounds). Set `use_dummy_constant = true`
to work around this: when `l` or `r` is nonzero it appends a variable pinned to
`1` (via a linear equality, so all subsolvers work), folds `l`, `r` into an extra
column of `L`, `R`, and uses `lbL = lbR = 0`. The extra variable is dropped from
the returned `x`.

Returns a NamedTuple:

    (; converged, x, y, objective, iterations, outer_iterations,
       subproblem_iterations, rho, status, message, return_value,
       solve_time_seconds)
"""
function solve_qpcc_with_lcqpow(data::NamedTuple;
                                x0=nothing,
                                print_level::Integer=PRINT_NONE,
                                use_dummy_constant::Bool=false,
                                kwargs...)
    q = data.q
    n = length(q)
    Qm = data.Q
    c0 = data.c0

    # Stack equalities then inequalities into LCQPow's two-sided linear block:
    #   J_eq x + b_eq == 0       ->  lbA = ubA = -b_eq
    #   J_ineq x + b_ineq >= 0   ->  lbA = -b_ineq,  ubA = +∞
    # Zero-row blocks collapse to empty arrays, which LCQPow reads as "no rows".
    J_eq   = data.J_eq
    b_eq   = data.b_eq
    J_ineq = data.J_ineq
    b_ineq = data.b_ineq
    A   = vcat(J_eq, J_ineq)
    lbA = vcat(-b_eq, -b_ineq)
    ubA = vcat(-b_eq, fill(INFTY, length(b_ineq)))

    Lm = data.L
    nComp = size(Lm, 1)
    Rm = data.R
    l = data.l
    r = data.r
    x0v = isnothing(x0) ? Float64[] : Vector{Float64}(x0)

    # Complementarity `0 <= Lx + l  ⊥  Rx + r >= 0` maps to `lbL = -l, lbR = -r`.
    # LCQPow can be unreliable with nonzero complementarity bounds, so when
    # `use_dummy_constant` is set and l/r are not all zero, append a variable
    # pinned to 1 and fold l, r into an extra column of L, R, leaving lbL = lbR
    # = 0. The dummy is pinned with a linear equality row (rather than a box
    # bound, which OSQP rejects). It stays out of the objective and is dropped
    # from the returned solution.
    if use_dummy_constant && (any(!iszero, l) || any(!iszero, r))
        Qs   = [Qm zeros(n); zeros(1, n + 1)]
        gs   = vcat(q, 0.0)
        Ls   = hcat(Lm, l)
        Rs   = hcat(Rm, r)
        lbLs, lbRs = zeros(nComp), zeros(nComp)
        pin  = reshape(vcat(zeros(n), 1.0), 1, n + 1)   # e_dummy' x == 1
        As   = size(A, 1) == 0 ? pin : vcat(hcat(A, zeros(size(A, 1))), pin)
        lbAs = vcat(lbA, 1.0)
        ubAs = vcat(ubA, 1.0)
        x0s  = isempty(x0v) ? Float64[] : vcat(x0v, 1.0)
    else
        Qs, gs, Ls, Rs, As = Qm, q, Lm, Rm, A
        lbAs, ubAs = lbA, ubA
        lbLs, lbRs = -l, -r
        x0s = x0v
    end

    options = _option_values((; print_level=print_level, kwargs...))

    res = _solve_qpcc_with_lcqpow(Qs, gs, Ls, Rs,
                                  lbLs, Float64[], lbRs, Float64[],
                                  As, lbAs, ubAs,
                                  Float64[], Float64[],
                                  x0s, Float64[],
                                  options)

    # Drop the dummy variable (if any) and score with the original objective.
    x = Vector{Float64}(primal_solution(res))[1:n]
    objective = 0.5 * dot(x, Qm * x) + dot(q, x) + c0

    # LCQPow has no message string for a successful return; supply a sensible one.
    msg = return_value(res) == 0 ? "Successful return" : _clean(message(res))

    return (
        converged             = converged(res),
        x                     = x,
        y                     = Vector{Float64}(dual_solution(res)),
        objective             = objective,
        iterations            = Int(iter_total(res)),
        outer_iterations      = Int(iter_outer(res)),
        subproblem_iterations = Int(subproblem_iter(res)),
        rho                   = rho_opt(res),
        status                = _clean(status_string(res)),
        message               = msg,
        return_value          = Int(return_value(res)),
        solve_time_seconds    = solve_time_seconds(res),
    )
end

end # module
