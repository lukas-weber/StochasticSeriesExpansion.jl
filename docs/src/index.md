# StochasticSeriesExpansion.jl

An implementation of the stochastic series expansion (SSE) quantum Monte Carlo algorithm for Julia.

## Features

* Run high-performance SSE simulations without being a QMC expert
* Supports anisotropic spin-``S`` quantum magnets out of the box
* Extensible to arbitrary models (but may suffer from sign problem)

## Getting started

Jump right into the Tutorial (TODO).

## Can I simulate my own model?

Yes, by implementing the [AbstractModel](@ref abstract_model) interface, which requires you the graph of bonds and arbitrary bond operators that describe your Hamiltonian. A key ingredient of the implementation is the “abstract loop” algorithm which will automatically figure out a set of worm updates (which is traditionally choosen by hand). How well the algorithm performs in practice depends on the model, in particular if there is a sign problem. So far it has been applied to

* [Heisenberg-type magnets with various anisotropic terms in the single-site ``S^z`` basis](https://doi.org/10.1038/s42005-023-01359-x)
* [Frustrated Heisenberg magnets in multi-site bases](https://doi.org/10.21468/scipostphys.12.2.054)
* [Magnets coupled to cavity photons](https://doi.org/10.1103/physrevb.100.054437)
