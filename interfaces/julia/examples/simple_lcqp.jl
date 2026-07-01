using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using LCQPow

# minimize    1/2 (x1-1)^2 + 1/2 (x2-1)^2
#           = 1/2 x'x - [1,1]'x + const
# subject to  0 <= x1  complements  x2 >= 0
#
# The complementary optimum is one of (1,0) or (0,1).
data = (
    Q = [1.0 0.0
         0.0 1.0],
    q = [-1.0, -1.0],
    c0 = 0.0,
    J_eq   = zeros(0, 2), b_eq   = zeros(0),   # no equality constraints
    J_ineq = zeros(0, 2), b_ineq = zeros(0),   # no inequality constraints
    L = [1.0 0.0], l = [0.0],                  # L*x + l = x1
    R = [0.0 1.0], r = [0.0],                  # R*x + r = x2
)

# Pick the QP subsolver via `qp_solver`:
#   LCQPow.QPOASES_DENSE   (default), LCQPow.QPOASES_SPARSE, LCQPow.OSQP_SPARSE
for (name, backend) in (("qpOASES", LCQPow.QPOASES_DENSE), ("OSQP", LCQPow.OSQP_SPARSE))
    result = LCQPow.solve_qpcc_with_lcqpow(data;
        qp_solver = backend,
        stationarity_tolerance = 1e-8,
        complementarity_tolerance = 1e-8,
    )

    println("\n=== $name ===")
    @show result.converged
    @show result.x
    @show result.objective
    @show result.iterations
    @show result.status
    @show result.solve_time_seconds
end
