# [OpstringEstimator interface](@id abstract_opstring_estimator)

The `AbstractOpstringEstimator` interface is used to implement most observable estimators in StochasticSeriesExpansion.

```@docs
AbstractOpstringEstimator
StochasticSeriesExpansion.init
StochasticSeriesExpansion.measure
StochasticSeriesExpansion.result
StochasticSeriesExpansion.register_evaluables
```

The reference implementation of this is the general purpose `MagnetizationEstimator` which works for all models that have
some kind of magnetization.
```@docs
MagnetizationEstimator
magnetization_state
magnetization_lattice_site_idx
magest_standard_prefix
```
