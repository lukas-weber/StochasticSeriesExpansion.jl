# Simulation data structure

In order to be generic, StochasticSeriesExpansion.jl internally uses data structures that
do not know a lot about physics but retain only the necessary information to run the SSE algorithm.
As a practicioner, you need to know about these structures either.

If you want to implement your own models or measurements, however, it is your job to translate physics
into these data structures, an `SSEData` object to be precise.

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
