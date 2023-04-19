using LoadLeveller

abstract type AbstractEstimator end

macro stub(func::Expr)
    return :(
        $func = error(
            "AbstractEstimator interface not implemented for Estimator of type $(typeof(est))",
        )
    )
end

@stub init(est::AbstractEstimator, state::AbstractVector{<:StateIdx})

@stub measure(
    est::AbstractEstimator,
    op::OperCode,
    state::AbstractVector{<:StateIdx},
    sse_data::SSEData,
)

@stub result(est::AbstractEstimator, ctx::MCContext)
