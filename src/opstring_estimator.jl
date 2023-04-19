using LoadLeveller

abstract type AbstractOpstringOpstringEstimator end

macro stub(func::Expr)
    return :(
        $func = error(
            "AbstractOpstringEstimator interface not implemented for Estimator of type $(typeof(est))",
        )
    )
end

@stub init(est::AbstractOpstringEstimator, state::AbstractOpstringVector{<:StateIdx})

@stub measure(
    est::AbstractOpstringEstimator,
    op::OperCode,
    state::AbstractOpstringVector{<:StateIdx},
    sse_data::SSEData,
)

@stub result(est::AbstractOpstringEstimator, ctx::MCContext)
