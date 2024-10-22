"""
The type used to define models that can be simulated by StochasticSeriesExpansion.jl.

Models are expected to have a constructor

    YourModel(parameters::AbstractDict{Symbol, <:Any})

that gets passed the [Carlo.jl](https://github.com/lukas-weber/Carlo.jl) task parameters that can be used to describe all parameters of your model, such as the lattice or the coupling strengths.

Apart from this, methods for the following functions should be implemented:

- [`generate_sse_data`](@ref)
- [`get_opstring_estimators`](@ref)
- [`leg_count`](@ref)
- [`normalization_site_count`](@ref)
"""
abstract type AbstractModel end

"""
    leg_count(model::Type{YourModel}) -> Integer

This function returns the *maximum* number of legs a bond operator can have in the model.

In the SSE algorithm, the Hamiltonian is decomposed into *bond operators* that (ideally) act on only a few sites. For example, the Heisenberg model consists of operators that act on two sites. In a diagrammatic picture, each site corresponds to one incoming and one outgoing leg. Therefore, the leg count in the Heisenberg model is four.

Operators with differing leg-counts are supported within the same model. In such cases, `leg_count` should return the maximum.
"""
function leg_count end

"""
    generate_sse_data(model::YourModel) -> SSEData

This function should translate the data saved in the `model` to construct the information needed by the abstract-loop SSE algorithm:

- a graph of bonds
- the bond Hamiltonians

From this information an [`SSEData`](@ref) instance can be constructed and returned.
"""
function generate_sse_data end

"""
    normalization_site_count(model::YourModel) -> Integer

The number of physical sites used for normalizing observables. It may differ from the SSE algorithmic site count sometimes, e.g. when the computational basis consists of multiple physical spins but we still want the energy to be measured per spin.
"""
function normalization_site_count end

"""
    get_opstring_estimators(model::YourModel) -> Vector{Type{<:AbstractOpstringEstimator}}

Returns an array of types of operator string estimators that are used to calculate most observables. Each of them should implement the [Operator string estimator](@ref abstract_opstring_estimator) interface.
"""
function get_opstring_estimators end
