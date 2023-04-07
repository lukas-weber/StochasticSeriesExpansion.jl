using LoadLeveller
using HDF5

using Random

Base.@kwdef mutable struct MC{Model <: AbstractModel} <: LoadLeveller.AbstractMC
    T::Float64 = 0.0

    avg_worm_length::Float64 = 1.0
    num_worms::Float64 = 0.0

    num_operators::Int64 = 0.0

    operators::Vector{OperCode}
    state::Vector{StateIdx}

    model::Model
    sse_data::SSEData

    vertices::Matrix{Int}
    v_first::Vector{Int}
    v_last::Vector{Int}
end

function MC{Model}(params::AbstractDict)
    model = Model(params)
    sse_data = generate_sse_data(model)
    return MC(
        v_first = zeros(Int, site_count(sse_data)),
        v_last = zeros(Int, site_count(sse_data)),
        T = params[:T]
    )
end

function LoadLeveller.init!(mc::MC, ctx::LoadLeveller.MCContext, params::AbstractDict)
    mc.state = [floor(StateIdx, rand(ctx.rng, site_count(mc.sse_data)) * get_site_data(mc.sse_data, i).dim) for i in 1:site_count(mc.sse_data)]

    init_opstring_cutoff = get(params, :init_opstring_cutoff, round(Int, site_count(mc.sse_data) * mc.T))
    mc.operators = [OperCode() for i in 1:init_opstring_cutoff]

    diagonal_warmup_sweeps = get(params, :diagonal_warmup_sweeps, 5)
    for i in 1:diagonal_warmup_sweeps
        diagonal_update(mc, ctx)
    end
end

function LoadLeveller.sweep!(mc::MC, ctx::LoadLeveller.MCContext)
    diagonal_update(mc, ctx)
    make_vertex_list(mc)
    worm_update()
end

function LoadLeveller.measure!(mc::MC, ctx::LoadLeveller.MCContext)
    sign = measure_sign(mc)

    measure!(ctx, :Sign, sign)
    measure!(ctx, :OperatorCount, mc.num_operators)
    measure!(ctx, :SignOperatorCount, sign*mc.num_operators)
    measure!(ctx, :SignOperatorCount2, sign*Float64(mc.num_operators)^2)
    measure!(ctx, :SignEnergy, -sign*mc.num_operators*mc.T-mc.sse_data.energy_offset / normalization_site_count(mc.model))
end
    