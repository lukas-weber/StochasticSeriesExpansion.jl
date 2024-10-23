# Interfaces

This page lists the different interfaces that need to be implemented to run StochasticSeriesExpansion.jl with custom models or estimators. Their purpose is to translate the information of a physical model into a representation that the SSE algorithm can understand.

## [Model interface](@id abstract_model)

The `AbstractModel` model interface describes a Hamiltonian that can be simulated with StochasticSeriesExpansion.jl. See [Available models](@ref available_models) to see the available example implementations. 

```@docs
AbstractModel
StochasticSeriesExpansion.generate_sse_data
StochasticSeriesExpansion.get_opstring_estimators
StochasticSeriesExpansion.leg_count
StochasticSeriesExpansion.normalization_site_count
```

## [Operator string estimator interface](@id abstract_opstring_estimator)

The `AbstractOpstringEstimator` interface is used to implement most observable estimators in StochasticSeriesExpansion.jl.

```@docs
AbstractOpstringEstimator
StochasticSeriesExpansion.init
StochasticSeriesExpansion.measure
StochasticSeriesExpansion.result
StochasticSeriesExpansion.register_evaluables
```

### MagnetizationEstimator
A useful general purpose implementation of the operator string estimator interface is the `MagnetizationEstimator` which works for all models that have
some kind of magnetization.
```@docs
MagnetizationEstimator
magnetization_state
magnetization_lattice_site_idx
magnetization_estimator_standard_prefix
```
