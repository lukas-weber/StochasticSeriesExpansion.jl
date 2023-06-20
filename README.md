# StochasticSeriesExpansion

This is a state-of-the-art implementation of the stochastic series expansion quantum Monte Carlo (QMC) algorithm [[1]](#1) with abstract loop updates [[2]](#2), which allow
the efficient simulation of many different bosonic lattice models in different computational bases using the same code.

StochasticSeriesExpansion is aimed both at
* **Nonspecialists** who want to produce unbiased QMC data for comparison against experiments or other methods,
and
* **Specialists** who want to extend it to study novel models and phenomena.

Out of the box, it can simulate common quantum magnet Hamiltonians, but it can be easily extended to arbitrary models by providing a bond Hamiltonian matrix.

It is built on the [LoadLeveller](https://github.com/lukas-weber/LoadLeveller.jl.git) framework.

## Getting started

To start, install LoadLeveller and StochasticSeriesExpansion in Julia
```julia
using Pkg
Pkg.add("LoadLeveller")
Pkg.add("StochasticSeriesExpansion")
```

Next, follow the [tutorial](example/TUTORIAL.md) to see how to simulate a Heisenberg model with single-ion anisotropies.

## References
<a id="1">[1]</a> A. W. Sandvik, Phys. Rev. B **59**, R14157(R) (1999)

<a id="2">[2]</a> L. Weber, A. Honecker, B. Normand, P. Corboz, F. Mila, S. Wessel, SciPost Phys. **12**, 054 (2022)
