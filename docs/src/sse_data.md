# Simulation data structure

In order to be generic, StochasticSeriesExpansion.jl internally uses data structures that
do not know a lot about physics but retain only the necessary information to run the SSE algorithm.

If you want to implement your own models or measurements, you have to translate physics
into into an `SSEData` object.

```@docs
SSEData
get_vertex_data
```

To construct an `SSEData` object you need to assemble the `VertexData` and `SSEBond` objects.

```@docs
SSEBond
VertexData
get_leg_state
```
