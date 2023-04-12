abstract type AbstractModel end

macro stub(func::Expr)
    return :(
        $func = error(
            "AbstractModel interface not implemented for Model type $(typeof(model))",
        )
    )
end

macro stubT(func::Expr)
    return :(
        $func = error("AbstractModel interface not implemented for Model type $(Model)")
    )
end

@stubT leg_count(::Type{<:AbstractModel})
@stub generate_sse_data(model::AbstractModel)::SSEData
