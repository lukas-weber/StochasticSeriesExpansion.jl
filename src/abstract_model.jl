abstract type AbstractModel end

@stubT leg_count(model::Type{<:AbstractModel})
@stub generate_sse_data(model::AbstractModel)::SSEData
@stub normalization_site_count(model::AbstractModel)
