"""
The type used to define Models that can be simulated by StochasticSeriesExpansion. See [`StochasticSeriesExpansion.Models.Magnet`](@ref) for a reference implementation.
"""
abstract type AbstractModel end

"""
    leg_count(model::Type{<:AbstractModel})

This function returns the *maximum* number of legs a bond operator can have in the model.

In the SSE algorithm, the Hamiltonian is decomposed into *bond operators* that (ideally) act on only a few sites. For example, the Heisenberg model consists of operators that act on two sites. In a diagrammatic picture, each site corresponds to one incoming and one outgoing leg. Therefore, the leg count in the Heisenberg model is four.
"""
function leg_count end

"""
    generate_sse_data(model::AbstractModel) -> SSEData

This function is responsible for generating an [`SSEData`](@ref) instance for `model`. It should use the data saved in the model to construct the information needed by the SSE algorithm, such as the bonds, their connection, and the Hamiltonians acting on each bond.
"""
function generate_sse_data end

"""
    normalization_site_count(model::AbstractModel) -> Int

The number of physical sites count used for normalizing observables. It may differ from the SSE algorithmic site count sometimes, e.g. when the computational basis consists of multiple sites but we still want e.g. our energy to be measured per spin.
"""
function normalization_site_count end

"""
    get_opstring_estimators(model::AbstractModel) -> Vector{Type{<:AbstractOpstringEstimator}}

Returns an array of types of operator string estimators that are used to calculate most observables. Each of them should implement the [`AbstractOpstringEstimator`](@ref abstract_opstring_estimator) interface.
"""
function get_opstring_estimators end
