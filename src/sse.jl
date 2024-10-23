using Carlo
using HDF5

using Random

mutable struct MC{Model<:AbstractModel,NSites} <: AbstractMC
    T::Float64
    opstring_estimators::Vector{Type}

    target_worm_length_fraction::Float64
    avg_worm_length::Float64
    num_worms::Float64

    num_operators::Int64

    operators::Vector{OperCode}
    state::Vector{StateIndex}

    model::Model
    sse_data::SSEData{NSites}

    vertex_list::VertexList{NSites}
end

function MC(params::AbstractDict)
    model = params[:model](params)
    sse_data = generate_sse_data(model)

    nsites = leg_count(typeof(model)) ÷ 2
    return MC{typeof(model),nsites}(
        params[:T],
        get_opstring_estimators(model),
        get(params, :target_worm_length_fraction, 2.0),
        1.0,
        5.0,
        0,
        OperCode[],
        StateIndex[],
        model,
        sse_data,
        VertexList{nsites}(length(sse_data.sites)),
    )
end

function Carlo.init!(mc::MC, ctx::MCContext, params::AbstractDict)
    mc.state = [rand(ctx.rng, StateIndex.(1:s.dim)) for s in mc.sse_data.sites]

    init_opstring_cutoff =
        get(params, :init_opstring_cutoff, round(Int, length(mc.sse_data.sites) * mc.T))
    mc.operators = [OperCode(Identity) for _ = 1:init_opstring_cutoff]

    diagonal_warmup_sweeps = get(params, :diagonal_warmup_sweeps, 5)
    for _ = 1:diagonal_warmup_sweeps
        diagonal_update(mc, ctx)
    end

    return nothing
end

function Carlo.sweep!(mc::MC, ctx::MCContext)
    diagonal_update(mc, ctx)
    make_vertex_list!(mc.vertex_list, mc.operators, mc.sse_data.bonds)
    worm_update(mc, ctx)

    return nothing
end

function Carlo.measure!(mc::MC, ctx::MCContext)
    sign = measure_sign(mc.operators, mc.sse_data)

    measure!(ctx, :Sign, float(sign))
    measure!(ctx, :OperatorCount, float(mc.num_operators))
    measure!(ctx, :SignOperatorCount, float(sign * mc.num_operators))
    measure!(ctx, :SignOperatorCount2, sign * float(mc.num_operators)^2)
    measure!(
        ctx,
        :SignEnergy,
        -sign * (mc.num_operators * mc.T + mc.sse_data.energy_offset) /
        normalization_site_count(mc.model),
    )

    measure_opstring!(mc, ctx, sign, mc.opstring_estimators...)

    return nothing
end

function Carlo.write_checkpoint(mc::MC, out::HDF5.Group)
    out["num_operators"] = mc.num_operators
    out["avg_worm_length"] = mc.avg_worm_length
    out["num_worms"] = mc.num_worms
    out["operators"] = mc.operators
    out["state"] = mc.state

    return nothing
end

function Carlo.read_checkpoint(mc::MC, in::HDF5.Group)
    mc.num_operators = in["num_operators"]
    mc.avg_worm_length = in["avg_worm_length"]
    mc.num_worms = in["num_worms"]
    mc.operators = in["operators"]
    mc.state = in["state"]

    return nothing
end

unsign(signobs, sign) = signobs ./ sign

function Carlo.register_evaluables(::Type{<:MC}, eval::Evaluator, params::AbstractDict)
    model = params[:model](params)

    for estimator in get_opstring_estimators(model)
        register_evaluables(estimator, eval, params)
    end

    #register_evaluables(model, eval, params)

    evaluate!(unsign, eval, :Energy, (:SignEnergy, :Sign))
    evaluate!(
        eval,
        :SpecificHeat,
        (:SignOperatorCount2, :SignOperatorCount, :Sign),
    ) do sn2, sn, s
        return (sn2 / s - sn * sn / s^2 - sn / s) / normalization_site_count(model)
    end

    return nothing
end


function diagonal_update(mc::MC{Model,NSites}, ctx::MCContext) where {Model,NSites}
    if mc.num_operators >= length(mc.operators) * 0.5
        if is_thermalized(ctx)
            @warn "spin array resized after thermalization"
        end
        old_length = length(mc.operators)
        resize!(mc.operators, ceil(Int64, old_length * 1.5 + 100))
        mc.operators[old_length+1:end] .= [OperCode(Identity)]
    end

    p_make_bond_raw = length(mc.sse_data.bonds) / mc.T
    p_remove_bond_raw = mc.T / length(mc.sse_data.bonds)

    for (iop, op) in enumerate(mc.operators)
        if isidentity(op)
            bond = rand(ctx.rng, 1:length(mc.sse_data.bonds))
            b = mc.sse_data.bonds[bond]

            # reversed because of join_idx convention
            dims = tuple((mc.sse_data.sites[site].dim for site in b.sites)...)
            idxs = tuple((mc.state[site] for site in b.sites)...)
            state_idx = join_idx(dims, idxs)

            vertex_data = get_vertex_data(mc.sse_data, bond)
            new_vert = get_diagonal_vertex(vertex_data, state_idx)
            weight = get_vertex_weight(vertex_data, new_vert)

            p_make_bond = p_make_bond_raw / (length(mc.operators) - mc.num_operators)

            if rand(ctx.rng) < p_make_bond * weight
                mc.operators[iop] = OperCode(bond, new_vert)
                mc.num_operators += 1
            end
        else
            bond = get_bond(op)
            vertex_data = get_vertex_data(mc.sse_data, bond)

            if isdiagonal(op)
                weight = get_vertex_weight(vertex_data, get_vertex(op))
                p_remove_bond =
                    (length(mc.operators) - mc.num_operators + 1) * p_remove_bond_raw
                if rand(ctx.rng) * weight < p_remove_bond
                    mc.operators[iop] = OperCode(Identity)
                    mc.num_operators -= 1
                end
            else
                b = mc.sse_data.bonds[bond]
                leg_state = get_leg_state(vertex_data, get_vertex(op))
                for s = 1:NSites
                    mc.state[b.sites[s]] = leg_state[NSites+s]
                end
            end
        end
    end
end

function worm_update(mc::MC, ctx::MCContext)
    total_worm_length = 1.0
    for _ = 1:ceil(Int, mc.num_worms)
        worm_length = worm_traverse!(mc, ctx)
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
        (l, p) = mc.vertex_list.v_first[i]
        if p < 0
            mc.state[i] = rand(ctx.rng, 1:mc.sse_data.sites[i].dim)
        else
            op = mc.operators[p]
            mc.state[i] =
                get_leg_state(get_vertex_data(mc.sse_data, get_bond(op)), get_vertex(op))[l]
        end
    end

    return nothing
end

function worm_traverse!(mc::MC{Model}, ctx::MCContext) where {Model}
    if mc.num_operators == 0
        return 0
    end

    p0 = 0
    l0 = 0

    while true
        p0 = rand(ctx.rng, 1:length(mc.operators))
        l0 = rand(ctx.rng, 1:leg_count(Model))
        if mc.vertex_list.vertices[l0, p0][1] > 0
            break
        end
    end

    op0 = mc.operators[p0]
    site0 = mc.sse_data.bonds[get_bond(op0)].sites[site_of_leg(l0, leg_count(Model) ÷ 2)]
    wormfunc0 = rand(ctx.rng, 1:worm_count(mc.sse_data.sites[site0].dim))

    return worm_traverse!(
        (l0, p0, wormfunc0),
        mc.operators,
        mc.vertex_list.vertices,
        mc.sse_data,
        ctx.rng,
    )
end

function worm_traverse!(
    start::NTuple{3,<:Integer},
    operators::AbstractVector{<:OperCode},
    vertices::AbstractArray,
    sse_data::SSEData{NSites},
    rng::AbstractRNG,
) where {NSites}
    (l0, p0, wormfunc0) = start
    (leg_in, p, wormfunc) = start

    worm_length = 1

    while true
        op = operators[p]
        bond = get_bond(op)
        (leg_out, wormfunc_out, new_vertex) = scatter(
            get_vertex_data(sse_data, bond),
            get_vertex(op),
            leg_in,
            wormfunc,
            rand(rng),
        )

        operators[p] = OperCode(bond, new_vertex)
        site_out = sse_data.sites[sse_data.bonds[bond].sites[site_of_leg(leg_out, NSites)]]

        if p == p0 && leg_out == l0 && wormfunc_out == worm_inverse(wormfunc0, site_out.dim)
            break
        end

        worm_length += 1

        wormfunc = wormfunc_out
        (leg_in, p) = vertices[leg_out, p]

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
            sign += get_sign(get_vertex_data(data, get_bond(op)), get_vertex(op)) < 0
        end
    end

    return Bool(sign & 1) ? -1.0 : 1.0
end

"""
    measure_opstring!

This generated function fuses multiple SSE estimators that have to loop over the operator string into a single loop.
"""
@generated function measure_opstring!(
    mc::MC{Model},
    ctx::MCContext,
    sign::AbstractFloat,
    estimator_types::Type{<:AbstractOpstringEstimator}...,
) where {Model}
    get_type(::Type{Type{T}}) where {T} = T
    inits = Expr(
        :tuple,
        (:(init($(get_type(type)), mc.model, mc.state)) for type in estimator_types)...,
    )
    measures = Expr(
        :block,
        (
            :(measure(estimators[$i], op, mc.state, mc.sse_data)) for
            i = 1:length(estimator_types)
        )...,
    )
    results = Expr(
        :block,
        (:(result(estimators[$i], ctx, mc.T, sign)) for i = 1:length(estimator_types))...,
    )

    return quote
        estimators = $inits
        n = zero(Int64)
        for op in mc.operators
            if isidentity(op)
                continue
            end

            if !isdiagonal(op)
                b = mc.sse_data.bonds[get_bond(op)]
                vd = get_vertex_data(mc.sse_data, get_bond(op))

                leg_state = get_leg_state(vd, get_vertex(op))

                for (i, s) in enumerate(b.sites)
                    mc.state[s] = leg_state[leg_count(Model)÷2+i]
                end
            end

            if (n < mc.num_operators)
                $measures
            end

            n += 1
        end

        $results
        return nothing
    end
end

function print_opstring(operators::AbstractVector{<:OperCode}, data::SSEData)
    for (p, op) in enumerate(operators)
        if !isidentity(op)
            println(
                "$(p): $(data.bonds[get_bond(op)].sites) - $(get_vertex_idx(get_vertex(op)))[$(join(Int.(get_leg_state(get_vertex_data(data, get_bond(op)), get_vertex(op))),","))]",
            )
        end
    end

    return nothing
end
