# StochasticSeriesExpansion.jl

An implementation of the stochastic series expansion (SSE) quantum Monte Carlo algorithm for Julia.

## Features

* Run high-performance SSE simulations without being a QMC expert
* Supports anisotropic spin-``S`` quantum Magnets out of the box
* Extensible to arbitrary models (but may suffer from sign problem)

## Getting started

Jump right into the Tutorial (TODO).

## Can I simulate my own model?

All you need to implement is the [AbstractModel](@ref abstract_model) interface. How well the algorithm performs in the end, depends on the model, in particular if there is a sign problem. Feel free to reach out (e.g. on Github) for assistance!
