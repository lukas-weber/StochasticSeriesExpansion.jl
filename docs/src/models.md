# [Predefined models](@id available_models)

This part of the documentation is about the predefined implementations of the [model interface](@ref abstract_model) that are included with StochasticSeriesExpansion.jl. At the moment there is one such implementation.

## MagnetModel

The `MagnetModel` implementation allows the simulation of Heisenberg-type magnets on arbitrary lattices.

```@docs
StochasticSeriesExpansion.MagnetModel
```

## ClusterModel

!!! warning
    This feature is experimental. It is tested and should work, but its API and parameters may change significantly in a future release.

```@docs
StochasticSeriesExpansion.ClusterModel
StochasticSeriesExpansion.ClusterBasis
StochasticSeriesExpansion.ClusterBases
```

## Common ingredients
### Lattice
This is a simple lattice structure that may be shared across different models. Practitioners only need to know how to create lattices from the task parameters.
```@docs
StochasticSeriesExpansion.Lattice
StochasticSeriesExpansion.UnitCell
StochasticSeriesExpansion.UCSite
StochasticSeriesExpansion.UCBond
StochasticSeriesExpansion.UnitCells
```
The following are internal pieces of the lattice structure that are useful for model implementers.
```@docs
StochasticSeriesExpansion.LatticeSite
StochasticSeriesExpansion.LatticeBond
```
### [Parameter maps](@id parameter_maps)
By default, model parameters like the exchange coupling `J` of the `MagnetModel` are equal for all bonds in the unit cell. Parameter maps are a way to override this behavior and assign different parameters to different unit cell positions.

The following example sets different values for the parameter `J` in the ``x`` and ``y`` direction.
```@example
using Carlo.JobTools
using StochasticSeriesExpansion

tm = TaskMaker()
tm.model = MagnetModel
tm.lattice = (;
    UnitCells.square,
)

tm.parameter_map = (;
    J = (:Jx, :Jy)
)

tm.Jx = 1.0
tm.Jy = 0.5
```
The length and order of the tuple `(:Jx, :Jy)` is tied to the way the bonds are specified in the unitcell structure `UnitCells.square`. In a similar way, one can change site-properties such as `h` or `S` by referring to the different sites in the unit cell.
```@example
using Carlo.JobTools
using StochasticSeriesExpansion

tm = TaskMaker()
tm.model = MagnetModel
tm.lattice = (;
    UnitCells.honeycomb,
)

tm.parameter_map = (;
    J = (:J1, :J2, :J3),
    S = (:Sa, :Sb),
)

tm.J1 = 1.0
tm.J2 = 0.5
tm.J3 = 1.5

tm.Sa = 1//2
tm.Sb = 1
```
The above example specifies a honeycomb lattice with spin-1/2 on one sublattice and spin-1 on the other. The three directions of nearest-neighbor exchange interactions have different strengths.
