using LoadLeveller

abstract type AbstractOpstringEstimator end

@stub init(est::AbstractOpstringEstimator, state::AbstractVector{<:StateIdx})

@stub measure(
    est::AbstractOpstringEstimator,
    op::OperCode,
    state::AbstractVector{<:StateIdx},
    sse_data::SSEData,
)

@stub result(
    est::AbstractOpstringEstimator,
    ctx::MCContext,
    T::AbstractFloat,
    Sign::AbstractFloat,
)
@stubT register_evaluables(
    est::Type{<:AbstractOpstringEstimator},
    eval::LoadLeveller.Evaluator,
)
