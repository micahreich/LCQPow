# LCQPow Julia interface

Minimal [CxxWrap](https://github.com/JuliaInterop/CxxWrap.jl) bindings that let
you call LCQPow's solver from Julia.

## Build

From the repository root, configure with the Julia interface enabled and build
the `lcqpow_julia` target (this also builds the qpOASES / OSQP dependencies):

```bash
cmake -S . -B build -DBUILD_JULIA_INTERFACE=ON
cmake --build build --target lcqpow_julia -j4
```

`JlCxx_DIR` is auto-detected from your Julia installation's `CxxWrap` package
(install it first with `julia -e 'using Pkg; Pkg.add("CxxWrap")'`). The compiled
library is written to `build/lib/liblcqpow_julia.*`, which is where the Julia
package looks for it by default. If you build elsewhere, point it there with a
`LocalPreferences.toml` entry `liblcqpow_julia_path`.

Then instantiate the Julia project:

```bash
julia --project=interfaces/julia -e 'using Pkg; Pkg.instantiate()'
```

## Usage

The solver takes QPCC data in CRISP's "Marble" form as a `NamedTuple`:

```
minimize    1/2 x'Q x + q'x + c0
subject to  J_eq x + b_eq == 0
            J_ineq x + b_ineq >= 0
            0 <= L x + l   complements   R x + r >= 0
```

```julia
using Pkg; Pkg.activate("interfaces/julia")
using LCQPow

# minimize  1/2 x'x - [1,1]x   s.t.   0 <= x1  complements  x2 >= 0
data = (
    Q = [1.0 0.0; 0.0 1.0],
    q = [-1.0, -1.0],
    c0 = 0.0,
    J_eq   = zeros(0, 2), b_eq   = zeros(0),   # no equality constraints
    J_ineq = zeros(0, 2), b_ineq = zeros(0),   # no inequality constraints
    L = [1.0 0.0], l = [0.0],                  # 0 <= x1
    R = [0.0 1.0], r = [0.0],                  #      complements  x2 >= 0
)

result = solve_qpcc_with_lcqpow(data;
    stationarity_tolerance = 1e-8,
    complementarity_tolerance = 1e-8,
)

result.converged   # true
result.x           # [1.0, 0.0]  (or [0.0, 1.0])
result.objective   # 1/2 x'Q x + q'x + c0
result.iterations
result.status      # e.g. "S-Stationary solution found"
```

Only `Q`, `q`, `L`, `l`, `R`, `r` are required; `c0`, `J_eq`, `b_eq`, `J_ineq`,
`b_ineq` are optional (omit them or pass zero-row blocks when there are no linear
constraints). An optional primal warm start is passed with the `x0` keyword.

`solve_qpcc_with_lcqpow` returns a `NamedTuple` with fields:

| field                   | meaning                                             |
| ----------------------- | --------------------------------------------------- |
| `converged`             | `true` iff LCQPow returned `SUCCESSFUL_RETURN`       |
| `x`                     | primal solution (length `n`)                        |
| `y`                     | dual solution (length depends on the subsolver)     |
| `objective`             | `1/2 x'Q x + q'x + c0` at the solution              |
| `iterations`            | total inner iterations                              |
| `outer_iterations`      | number of penalty updates                           |
| `subproblem_iterations` | total QP subsolver iterations                       |
| `rho`                   | penalty parameter at the solution                   |
| `status`                | stationarity type as a string                       |
| `message`               | human-readable return message                       |
| `return_value`          | raw LCQPow `ReturnValue` code (`0` = success)        |
| `solve_time_seconds`    | wall-clock time spent in `runSolver`                |

Internally this maps onto LCQPow's native form: equalities become
`lbA = ubA = -b_eq`, inequalities become `lbA = -b_ineq, ubA = +∞`, and the
complementarity maps to `lbL = -l`, `lbR = -r` (LCQPow enforces
`(Lx - lbL) ⊥ (Rx - lbR)`).

### Choosing the QP subsolver

`qp_solver` selects the backend LCQPow uses for its QP subproblems:

```julia
solve_qpcc_with_lcqpow(data; qp_solver = LCQPow.QPOASES_DENSE)   # default
solve_qpcc_with_lcqpow(data; qp_solver = LCQPow.QPOASES_SPARSE)
solve_qpcc_with_lcqpow(data; qp_solver = LCQPow.OSQP_SPARSE)     # OSQP
```

A sparse `qp_solver` (`QPOASES_SPARSE` or `OSQP_SPARSE`) automatically switches
the loaded problem to sparse mode.

### Nonzero complementarity constants

LCQPow can be unreliable when the complementarity constants `l`, `r` are nonzero
(they become nonzero `lbL`/`lbR` bounds). Pass `use_dummy_constant = true` to
apply the standard workaround — append a variable pinned to `1`, fold `l`, `r`
into an extra column of `L`, `R`, and use `lbL = lbR = 0`:

```julia
solve_qpcc_with_lcqpow(data; use_dummy_constant = true)
```

The trick only kicks in when `l` or `r` is actually nonzero (otherwise it's a
no-op), the dummy variable is dropped from the returned `x`, and it works with
every subsolver (it pins the dummy via a linear equality, not a box bound).

### Options

Other solver hyperparameters are keyword arguments (see `LCQPow.OPTION_NAMES`),
e.g. `max_iterations`, `initial_penalty_parameter`, `penalty_update_factor`, and
`print_level = LCQPow.PRINT_NONE | PRINT_OUTER | PRINT_INNER` (default
`PRINT_NONE`).

See `examples/simple_lcqp.jl` for a runnable example.
