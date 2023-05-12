using LoadLeveller

abstract type AbstractOpstringOpstringEstimator end

@stub init(est::AbstractOpstringEstimator, state::AbstractOpstringVector{<:StateIdx})

@stub measure(
    est::AbstractOpstringEstimator,
    op::OperCode,
    state::AbstractOpstringVector{<:StateIdx},
    sse_data::SSEData,
)

@stub result(est::AbstractOpstringEstimator, ctx::MCContext)
