using LoadLeveller
using HDF5

using Random

Base.@kwdef mutable struct MC{Model<:AbstractModel} <: LoadLeveller.AbstractMC
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

function MC{Model}(params::AbstractDict) where {Model}
    model = Model(params)
    sse_data = generate_sse_data(model)
    return MC(
        v_first = zeros(Int, site_count(sse_data)),
        v_last = zeros(Int, site_count(sse_data)),
        T = params[:T],
    )
end

function LoadLeveller.init!(mc::MC, ctx::LoadLeveller.MCContext, params::AbstractDict)
    mc.state = [
        floor(
            StateIdx,
            rand(ctx.rng, site_count(mc.sse_data)) * get_site_data(mc.sse_data, i).dim,
        ) for i = 1:site_count(mc.sse_data)
    ]

    init_opstring_cutoff =
        get(params, :init_opstring_cutoff, round(Int, site_count(mc.sse_data) * mc.T))
    mc.operators = [make_identity() for i = 1:init_opstring_cutoff]

    diagonal_warmup_sweeps = get(params, :diagonal_warmup_sweeps, 5)
    for i = 1:diagonal_warmup_sweeps
        diagonal_update(mc, ctx)
    end

    return nothing
end

function LoadLeveller.sweep!(mc::MC, ctx::LoadLeveller.MCContext)
    diagonal_update(mc, ctx)
    make_vertex_list(mc)
    worm_update()

    return nothing
end

function LoadLeveller.measure!(mc::MC, ctx::LoadLeveller.MCContext)
    sign = measure_sign(mc)

    measure!(ctx, :Sign, sign)
    measure!(ctx, :OperatorCount, mc.num_operators)
    measure!(ctx, :SignOperatorCount, sign * mc.num_operators)
    measure!(ctx, :SignOperatorCount2, sign * Float64(mc.num_operators)^2)
    measure!(
        ctx,
        :SignEnergy,
        -sign * mc.num_operators * mc.T -
        mc.sse_data.energy_offset / normalization_site_count(mc.model),
    )

    return nothing
end

function LoadLeveller.write_checkpoint(mc::MC, out::HDF5.Group)
    out["num_operators"] = mc.num_operators
    out["avg_worm_length"] = mc.avg_worm_length
    out["num_worms"] = mc.num_worms
    out["operators"] = mc.operators
    out["state"] = mc.state

    return nothing
end

function LoadLeveller.read_checkpoint(mc::MC, in::HDF5.Group)
    mc.num_operators = in["num_operators"]
    mc.avg_worm_length = in["avg_worm_length"]
    mc.num_worms = in["num_worms"]
    mc.operators = in["operators"]
    mc.state = in["state"]

    return nothing
end

function LoadLeveller.register_evaluables(
    ::Type{MC{Model}},
    eval::LoadLeveller.Evaluator,
    params::AbstractDict,
) where {Model}

    model = Model(params)
    register_evaluables(model, eval, params)

    evaluate!(eval, :Energy, [:SignEnergy, :Sign]) do energy, sign
        return energy / sign
    end

    evaluate!(
        eval,
        :SpecificHeat,
        [:SignOperatorCount2, :SignOperatorCount, :Sign],
    ) do sn2, sn, s
        return (sn2 / s - sn * sn / s^2 - sn / s) / normalization_site_count(model)
    end

    return nothing
end


function diagonal_update(mc::MC{Model}, ctx::LoadLeveller.MCContext) where {Model}
    if mc.num_operators >= length(mc.operators * 0.5)
        if is_thermalized(ctx)
            @warn "spin array resized after thermalization"
        end
        old_length = length(mc.operators)
        resize!(mc.operators, ceil(Int64, old_length * 1.5 + 100))
        mc.operators[old_length+1:end] .= OperCode(Identity)
    end

    p_make_bond_raw = bond_count(mc.sse_data) / mc.T
    p_remove_bond_raw = mc.T / bond_count(mc.sse_data)

    tmpstate = copy(mc.state)

    for (iop, op) in enumerate(mc.operators)
        if isidentity(op)
            bond = floor(Int, rand(ctx.rng) * bond_count(mc.sse_data))
            b = get_bond_data(mc.sse_data, bond)

            # TODO: define a function for this?
            state_idx = 0
            for s = 1:leg_count(Model)÷2
                state_idx *= get_site_data(mc.sse_data, b[s]).dim
                state_idx += mc.state[b[s]]
            end

            vertex_data = get_vertex_data(mc.sse_data, bond)
            new_vert = diagonal_vertex(vertex_data, state_idx)
            weight = vertex_weight(vertex_data, new_vert)

            p_make_bond = p_make_bond_raw / (length(mc.operators - mc.num_operators))

            if rand(ctx.rng) < p_make_bond * weight
                mc.operators[iop] = OperCode(bond, new_vert)
                mc.num_operators += 1
            end
        else
            bond = get_bond(op)
            vertex_data = get_vertex_data(mc.sse_data, bond)

            if isdiagonal(op)
                weight = vertex_weight(vertex_data, get_vertex(op))
                p_remove_bond =
                    (length(mc.operators) - mc.num_operators + 1) * p_remove_bond_raw
                if rand(ctx.rng) * weight < p_remove_bond
                    mc.operators[iop] = OperCode(Identity)
                    mc.num_operators -= 1
                end
            else
                b = get_bond_data(mc.sse_data, bond)
                leg_state = get_legstate(vertex_data, get_vertex(op))
                for s = 1:leg_count(Model)÷2
                    mc.state[b[s]] = leg_state[leg_count(Model)÷2+s]
                end
            end
        end
    end

    @assert tmpstate == mc.state
end
