# LCQPow

A solver for Quadratic Programs with Linear Complementarity Constraints.

This README gives the shortest path to building the C++ core and the Julia
bindings. For the full documentation (MATLAB/Python interfaces, algorithm
details, Docker, etc.) see [README.old.md](README.old.md).

## Prerequisites

- A C++ compiler and **CMake >= 3.13**
- **Julia >= 1.10** (only needed for the Julia bindings)

Install the `CxxWrap` Julia package once — CMake auto-detects its C++ library
(`JlCxx`) from this:

```bash
julia -e 'using Pkg; Pkg.add("CxxWrap")'
```

## 1. Clone

Clone the repo and pull in the bundled dependencies (qpOASES, OSQP, googletest,
pybind11) as submodules:

```bash
git clone https://github.com/micahreich/LCQPow.git
cd LCQPow
git submodule update --init --recursive
```

## 2. Build the C++ core + Julia bindings

From the repository root:

```bash
cmake -S . -B build -DBUILD_JULIA_INTERFACE=ON
cmake --build build -j4
```

This builds the LCQPow library (and its qpOASES / OSQP dependencies) plus the
`lcqpow_julia` bindings. The compiled bindings land at
`build/lib/liblcqpow_julia.*`, which is where the Julia package looks by default.

> If you don't need Julia, just drop `-DBUILD_JULIA_INTERFACE=ON`. Examples and
> unit tests build by default; add `-DBUILD_EXAMPLES=OFF -DUNIT_TESTS=OFF` to
> skip them.

## 3. Set up the Julia package

Instantiate the Julia project so it picks up its dependencies:

```bash
julia --project=interfaces/julia -e 'using Pkg; Pkg.instantiate()'
```

## 4. Run it

```bash
julia --project=interfaces/julia interfaces/julia/examples/simple_lcqp.jl
```

Or from the REPL:

```julia
using Pkg; Pkg.activate("interfaces/julia")
using LCQPow

# minimize  1/2 x'x - [1,1]x   s.t.   0 <= x1  complements  x2 >= 0
data = (
    Q = [1.0 0.0; 0.0 1.0],
    q = [-1.0, -1.0],
    c0 = 0.0,
    J_eq   = zeros(0, 2), b_eq   = zeros(0),
    J_ineq = zeros(0, 2), b_ineq = zeros(0),
    L = [1.0 0.0], l = [0.0],
    R = [0.0 1.0], r = [0.0],
)

result = solve_qpcc_with_lcqpow(data)
result.x           # [1.0, 0.0]  (or [0.0, 1.0])
result.objective
result.converged
```

See [interfaces/julia/README.md](interfaces/julia/README.md) for the full Julia
API (problem form, options, QP-subsolver choice, warm starts).

## License

GNU Lesser General Public License (v2.1). See [LICENSE](LICENSE).
