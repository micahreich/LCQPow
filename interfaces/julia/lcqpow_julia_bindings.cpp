/*
 *  Minimal Julia (CxxWrap) bindings for LCQPow.
 *
 *  Exposes a single entry point, `_solve_qpcc_with_lcqpow`, that accepts a dense
 *  LCQP / QPCC in LCQPow's native formulation and returns a wrapped result
 *  object with named getters (converged, primal_solution, iter_total, ...).
 *  The thin Julia wrapper in `src/LCQPow.jl` repackages this into a NamedTuple.
 *
 *  LCQP formulation solved:
 *      minimize    1/2 x' Q x + g' x
 *      subject to  lbA <= A x <= ubA          (linear constraints, optional)
 *                  lb  <=  x  <= ub            (box constraints, optional)
 *                  lbL <= L x <= ubL   complements   lbR <= R x <= ubR
 */

#include <jlcxx/jlcxx.hpp>

#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "LCQProblem.hpp"
#include "MessageHandler.hpp"

using namespace LCQPow;

namespace {

/** Result handed back to Julia. Exposed through per-field getters (registered in
 *  the module below) so the Julia side can assemble a clean NamedTuple rather
 *  than unpacking a flat, index-addressed array. */
struct LCQPowResult {
    bool converged = false;
    int return_value = 0;
    std::string message;
    int status = 0;
    std::string status_string;
    int iter_total = 0;
    int iter_outer = 0;
    int subproblem_iter = 0;
    double rho_opt = 0.0;
    double solve_time_seconds = 0.0;
    std::vector<double> x;
    std::vector<double> y;
};

/** Convert a Julia (column-major) matrix into a row-major buffer, the layout
 *  expected by LCQPow's dense loader. Returns an empty vector for empty input. */
std::vector<double> to_row_major(jlcxx::ArrayRef<double, 2> arr) {
    const long rows = jl_array_dim(arr.wrapped(), 0);
    const long cols = jl_array_dim(arr.wrapped(), 1);
    std::vector<double> row_major(static_cast<size_t>(rows) * static_cast<size_t>(cols));
    const double* col = arr.data();
    for (long i = 0; i < rows; ++i)
        for (long j = 0; j < cols; ++j)
            row_major[static_cast<size_t>(i) * cols + j] = col[static_cast<size_t>(j) * rows + i];
    return row_major;
}

/** Pointer to the data, or nullptr when empty (LCQPow treats null as "use default"). */
const double* ptr_or_null(const std::vector<double>& v) { return v.empty() ? nullptr : v.data(); }
const double* ptr_or_null(jlcxx::ArrayRef<double, 1> v) { return v.size() == 0 ? nullptr : v.data(); }

/** Copy a std::vector into a freshly allocated Julia array. */
jlcxx::Array<double> to_julia(const std::vector<double>& v) {
    jlcxx::Array<double> out;
    for (double value : v)
        out.push_back(value);
    return out;
}

/** Option slots. Order must match `OPTION_NAMES` in the Julia wrapper. A NaN in
 *  a slot means "leave at the LCQPow default"; any other value is applied. */
enum OptionIndex {
    OPT_STATIONARITY_TOLERANCE = 0,
    OPT_COMPLEMENTARITY_TOLERANCE,
    OPT_INITIAL_PENALTY_PARAMETER,
    OPT_PENALTY_UPDATE_FACTOR,
    OPT_SOLVE_ZERO_PENALTY_FIRST,
    OPT_PERTURB_STEP,
    OPT_MAX_ITERATIONS,
    OPT_MAX_PENALTY_PARAMETER,
    OPT_N_DYNAMIC_PENALTY,
    OPT_ETA_DYNAMIC_PENALTY,
    OPT_PRINT_LEVEL,
    OPT_STORE_STEPS,
    OPT_QP_SOLVER,
    NUM_OPTIONS
};

void apply_options(Options& options, jlcxx::ArrayRef<double, 1> v) {
    if (v.size() != NUM_OPTIONS)
        throw std::runtime_error("LCQPow option vector has unexpected length");

    auto is_set = [&](int i) { return !std::isnan(v[i]); };
    auto as_int = [&](int i) { return static_cast<int>(std::lround(v[i])); };

    if (is_set(OPT_STATIONARITY_TOLERANCE))    options.setStationarityTolerance(v[OPT_STATIONARITY_TOLERANCE]);
    if (is_set(OPT_COMPLEMENTARITY_TOLERANCE)) options.setComplementarityTolerance(v[OPT_COMPLEMENTARITY_TOLERANCE]);
    if (is_set(OPT_INITIAL_PENALTY_PARAMETER)) options.setInitialPenaltyParameter(v[OPT_INITIAL_PENALTY_PARAMETER]);
    if (is_set(OPT_PENALTY_UPDATE_FACTOR))     options.setPenaltyUpdateFactor(v[OPT_PENALTY_UPDATE_FACTOR]);
    if (is_set(OPT_SOLVE_ZERO_PENALTY_FIRST))  options.setSolveZeroPenaltyFirst(v[OPT_SOLVE_ZERO_PENALTY_FIRST] != 0.0);
    if (is_set(OPT_PERTURB_STEP))              options.setPerturbStep(v[OPT_PERTURB_STEP] != 0.0);
    if (is_set(OPT_MAX_ITERATIONS))            options.setMaxIterations(as_int(OPT_MAX_ITERATIONS));
    if (is_set(OPT_MAX_PENALTY_PARAMETER))     options.setMaxPenaltyParameter(v[OPT_MAX_PENALTY_PARAMETER]);
    if (is_set(OPT_N_DYNAMIC_PENALTY))         options.setNDynamicPenalty(as_int(OPT_N_DYNAMIC_PENALTY));
    if (is_set(OPT_ETA_DYNAMIC_PENALTY))       options.setEtaDynamicPenalty(v[OPT_ETA_DYNAMIC_PENALTY]);
    if (is_set(OPT_PRINT_LEVEL))               options.setPrintLevel(as_int(OPT_PRINT_LEVEL));
    if (is_set(OPT_STORE_STEPS))               options.setStoreSteps(v[OPT_STORE_STEPS] != 0.0);
    if (is_set(OPT_QP_SOLVER))                 options.setQPSolver(as_int(OPT_QP_SOLVER));
}

LCQPowResult solve_qpcc_with_lcqpow(jlcxx::ArrayRef<double, 2> Q,
                                    jlcxx::ArrayRef<double, 1> g,
                                    jlcxx::ArrayRef<double, 2> L,
                                    jlcxx::ArrayRef<double, 2> R,
                                    jlcxx::ArrayRef<double, 1> lbL,
                                    jlcxx::ArrayRef<double, 1> ubL,
                                    jlcxx::ArrayRef<double, 1> lbR,
                                    jlcxx::ArrayRef<double, 1> ubR,
                                    jlcxx::ArrayRef<double, 2> A,
                                    jlcxx::ArrayRef<double, 1> lbA,
                                    jlcxx::ArrayRef<double, 1> ubA,
                                    jlcxx::ArrayRef<double, 1> lb,
                                    jlcxx::ArrayRef<double, 1> ub,
                                    jlcxx::ArrayRef<double, 1> x0,
                                    jlcxx::ArrayRef<double, 1> y0,
                                    jlcxx::ArrayRef<double, 1> option_values) {
    const int nV = static_cast<int>(g.size());
    const int nComp = static_cast<int>(jl_array_dim(L.wrapped(), 0));
    const int nC = static_cast<int>(jl_array_dim(A.wrapped(), 0));

    // LCQPow's dense loader expects row-major storage; Julia is column-major.
    const std::vector<double> Qd = to_row_major(Q);
    const std::vector<double> Ld = to_row_major(L);
    const std::vector<double> Rd = to_row_major(R);
    const std::vector<double> Ad = to_row_major(A);

    LCQProblem lcqp(nV, nC, nComp);

    Options options;
    apply_options(options, option_values);
    lcqp.setOptions(options);

    LCQPowResult result;

    ReturnValue ret = lcqp.loadLCQP(ptr_or_null(Qd), ptr_or_null(g),
                                    ptr_or_null(Ld), ptr_or_null(Rd),
                                    ptr_or_null(lbL), ptr_or_null(ubL),
                                    ptr_or_null(lbR), ptr_or_null(ubR),
                                    ptr_or_null(Ad), ptr_or_null(lbA), ptr_or_null(ubA),
                                    ptr_or_null(lb), ptr_or_null(ub),
                                    ptr_or_null(x0), ptr_or_null(y0));

    // A sparse subsolver needs the internally stored data converted to sparse.
    if (ret == SUCCESSFUL_RETURN && options.getQPSolver() >= QPOASES_SPARSE)
        ret = lcqp.switchToSparseMode();

    if (ret == SUCCESSFUL_RETURN) {
        const auto t0 = std::chrono::steady_clock::now();
        ret = lcqp.runSolver();
        const auto t1 = std::chrono::steady_clock::now();
        result.solve_time_seconds = std::chrono::duration<double>(t1 - t0).count();

        // A primal iterate exists whenever the load succeeded, even if the solve
        // itself stopped early (e.g. max iterations), so this is always safe here.
        result.x.assign(static_cast<size_t>(nV), 0.0);
        lcqp.getPrimalSolution(result.x.data());

        const int nDuals = lcqp.getNumberOfDuals();
        if (nDuals > 0) {
            result.y.assign(static_cast<size_t>(nDuals), 0.0);
            lcqp.getDualSolution(result.y.data());
        }

        OutputStatistics stats;
        lcqp.getOutputStatistics(stats);
        result.iter_total = stats.getIterTotal();
        result.iter_outer = stats.getIterOuter();
        result.subproblem_iter = stats.getSubproblemIter();
        result.rho_opt = stats.getRhoOpt();
        result.status = static_cast<int>(stats.getSolutionStatus());
        result.status_string = MessageHandler::SolutionString(stats.getSolutionStatus());
    }

    result.converged = (ret == SUCCESSFUL_RETURN);
    result.return_value = static_cast<int>(ret);
    result.message = MessageHandler::MessageString(ret);
    return result;
}

}  // namespace

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.add_type<LCQPowResult>("LCQPowResultCxx")
        .method("converged", [](const LCQPowResult& r) { return r.converged; })
        .method("return_value", [](const LCQPowResult& r) { return r.return_value; })
        .method("message", [](const LCQPowResult& r) { return r.message; })
        .method("status", [](const LCQPowResult& r) { return r.status; })
        .method("status_string", [](const LCQPowResult& r) { return r.status_string; })
        .method("iter_total", [](const LCQPowResult& r) { return r.iter_total; })
        .method("iter_outer", [](const LCQPowResult& r) { return r.iter_outer; })
        .method("subproblem_iter", [](const LCQPowResult& r) { return r.subproblem_iter; })
        .method("rho_opt", [](const LCQPowResult& r) { return r.rho_opt; })
        .method("solve_time_seconds", [](const LCQPowResult& r) { return r.solve_time_seconds; })
        .method("primal_solution", [](const LCQPowResult& r) { return to_julia(r.x); })
        .method("dual_solution", [](const LCQPowResult& r) { return to_julia(r.y); });

    mod.method("_solve_qpcc_with_lcqpow", &solve_qpcc_with_lcqpow);
}
