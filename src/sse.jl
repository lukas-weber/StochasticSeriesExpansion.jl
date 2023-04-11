using LoadLeveller
using HDF5

using Random

Base.@kwdef mutable struct MC{Model<:AbstractModel} <: LoadLeveller.AbstractMC
    T::Float64 = 0.0

    target_worm_length_fraction::Float64
    avg_worm_length::Float64 = 1.0
    num_worms::Float64 = 0.0

    num_operators::Int64 = 0.0

    operators::Vector{OperCode}
    state::Vector{StateIdx}

    model::Model
    sse_data::SSEData{Model}

    vertex_list::VertexList
end

function MC{Model}(params::AbstractDict) where {Model}
    model = Model(params)
    sse_data = generate_sse_data(model)
    return MC{Model}(
        vertex_list = VertexList(site_count(sse_data)),
        v_last = zeros(Int, site_count(sse_data)),
        T = params[:T],
        target_worm_length_fraction = get(params, :target_worm_length_fraction, 2.0),
        model = model,
        sse_data = sse_data,
    )
end

function LoadLeveller.init!(mc::MC, ctx::LoadLeveller.MCContext, params::AbstractDict)
    mc.state = [
        floor(StateIdx, rand(ctx.rng, site_count(mc.sse_data)) * mc.sse_data.sites[i].dim) for i = 1:site_count(mc.sse_data)
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
    make_vertex_list!(mc.vertex_list, mc.operators, mc.sse_data)
    worm_update()

    return nothing
end

function LoadLeveller.measure!(mc::MC, ctx::LoadLeveller.MCContext)
    sign = measure_sign(mc.operators, mc.sse_data)

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
            b = mc.sse_data.bonds[bond]

            # TODO: define a function for this?
            state_idx = 0
            for s = 1:leg_count(Model)÷2
                state_idx *= mc.sse_data.sites[b.sites[s]].dim
                state_idx += mc.state[b.sites[s]]
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
                b = mc.sse_data.bonds[bond]
                leg_state = get_legstate(vertex_data, get_vertex(op))
                for s = 1:leg_count(Model)÷2
                    mc.state[b.sites[s]] = leg_state[leg_count(Model)÷2+s]
                end
            end
        end
    end

    @assert tmpstate == mc.state
end

function worm_update(mc::MC, ctx::LoadLeveller.MCContext)
    total_worm_length = 1.0
    for i = 1:ceil(Int, mc.num_worms)
        worm_length = worm_traverse(mc, ctx)
        total_worm_length += worm_length
    end

    if is_thermalized(ctx) && mc.num_operators != 0
        measure!(ctx, :WormLengthFraction, total_worm_length / mc.num_operators)
    end

    avg_worm_length = total_worm_length / ceil(mc.num_worms)
    if !is_thermalized(ctx)
        mc.avg_worm_length += 0.01 * (avg_worm_length - mc.avg_worm_length)
        target_worms =
            mc.target_worm_length_fraction * mc.num_operators / mc.avg_worm_length

        mc.num_worms +=
            0.01 * (target_worms - mc.num_worms) + tanh(target_worms - mc.num_worms)
        mc.num_worms = clamp(mc.num_worms, 1.0, 1.0 + mc.num_operators / 2.0)
    end

    for i in eachindex(mc.state)
        (l, p) = mc.vertex_list.v_first[:, i]
        if p < 0
            mc.state[i] = rand(ctx.rng, 1:mc.sse_data.sites[i].dim)
        else
            op = mc.operators[p]
            mc.state[i] = get_vertex_data(mc.sse_data, get_bond(op)).leg_state[l]
        end
    end

    return nothing
end

function worm_traverse(mc::MC{Model}, ctx::LoadLeveller.MCContext) where {Model}
    if mc.num_operators == 0
        return 0
    end

    worm_length = 1

    p0 = rand(ctx.rng, 1:length(mc.operators))
    l0 = rand(ctx.rng, 1:leg_count(Model))

    op0 = mc.operators[p0]
    site0 = mc.sse_data.bonds[get_bond(op0)].sites[l0%(leg_count(Model)÷2)]
    wormfunc0 = rand(ctx.rng, 1:worm_count(mc.sse_data.sites[site0]))

    p = p0
    leg_in = l0
    wormfunc = wormfunc0

    while true
        op = mc.operators[p0]
        bond = get_bond(op)
        (leg_out, wormfunc_out, new_vertex) = scatter(
            get_vertex_data(
                mc.sse_data,
                bond,
                get_vertex(op),
                leg_in,
                wormfunc,
                rand(ctx.rng),
            ),
        )

        mc.operators[p] = OperCode(bond, new_vertex)
        site_out =
            mc.sse_data.sites[mc.sse_data.bonds[bond].sites[leg_out%(leg_count(Model)÷2)]]

        if p == p0 && leg_out == l0 && wormfunc_out == worm_inverse(wormfunc0, site_out)
            break
        end

        worm_length += 1

        wormfunc = wormfunc_out
        (leg_in, p) = mc.vertex_list.vertices[:, leg_out, p]

        if p == p0 && leg_in == l0 && wormfunc == wormfunc0
            break
        end
    end

    return worm_length
end

function measure_sign(operators::AbstractVector{<:OperCode}, data::SSEData)
    sign = 0
    for op in operators
        if !isidentity(op)
            sign += get_sign(get_vertex_data(data, get_bond(op)), op.vertex()) < 0
        end
    end

    return (sign & 1) ? -1 : 1
end
