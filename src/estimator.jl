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

function opstring_measurement(
    mc::MC{Model},
    ctx::MCContext,
    estimators::AbstractEstimator...,
) where {Model}
    for estimator in estimators
        init!(estimator, mc.state)
    end

    n = zero(Int64)

    for op in mc.operators
        if isidentity(op)
            continue
        end

        if !isdiagonal(op)
            b = mc.sse_data.bonds[get_bond(op)]
            vd = get_vertex_data(get_bond(op))

            leg_state = @view vd.leg_states[:, get_vertex(op)]

            mc.state[b.sites...] .= leg_state[leg_count(Model)÷2:end]
        end

        if (n < mc.num_operators)
            for estimator in estimators
                measure(estimator, op, mc.state, mc.sse_data)
            end
        end

        n += 1
    end

    for estimator in estimators
        result(estimator, ctx)
    end
end
