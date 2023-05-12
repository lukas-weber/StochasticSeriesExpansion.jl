using LoadLeveller

abstract type AbstractOpstringEstimator end

@stub init(est::AbstractOpstringEstimator, state::AbstractVector{<:StateIdx})

@stub measure(
    est::AbstractOpstringEstimator,
    op::OperCode,
    state::AbstractVector{<:StateIdx},
    sse_data::SSEData,
)

@stub result(est::AbstractOpstringEstimator, ctx::MCContext)
