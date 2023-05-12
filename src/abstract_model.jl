abstract type AbstractModel end

@stubT leg_count(model::Type{<:AbstractModel})
@stub generate_sse_data(model::AbstractModel)
@stub normalization_site_count(model::AbstractModel)
@stub get_opstring_estimators(model::AbstractModel, params::AbstractDict)
