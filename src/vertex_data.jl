using HiGHS
using LinearAlgebra
import MathOptInterface as MOI

Base.@kwdef struct Transition
    offset::Int = -1
    length::Int = 0
end

isinvalid(t::Transition) = t.offset < 0

struct VertexData{NSites}
    energy_offset::Float64
    dims::NTuple{NSites,Int}

    diagonal_vertices::Vector{VertexCode}
    signs::Vector{Int8}
    weights::Vector{Float64}

    transitions::Array{Transition,3} # [leg_in, worm_in, vertex]

    transition_cumprobs::Vector{Float64}
    transition_targets::Vector{VertexCode}
    transition_step_outs::Vector{Tuple{Int,Int}} # [leg, worm]

    leg_states::Matrix{StateIdx} # [leg, vertex]
end

function VertexData(
    dims::NTuple{NSites,<:Integer},
    bond_hamiltonian::AbstractMatrix;
    energy_offset_factor::AbstractFloat = 0.25,
) where {NSites}
    energy_offset = calc_energy_offset(bond_hamiltonian, energy_offset_factor)
    total_dim = prod(dims)
    max_worm_count = maximum(worm_count, dims)

    # TODO make variable
    tolerance = 1e-7
    lp_tolerance = 1e-10

    (diagonal_vertices, weights, leg_states, signs) =
        construct_vertices(dims, bond_hamiltonian, energy_offset, tolerance)

    @debug begin
        @assert NSites >= 1
        @assert total_dim == size(bond_hamiltonian, 1) == size(bond_hamiltonian, 2)
        @assert size(leg_states) == (2 * NSites, length(weights))
    end

    (transitions, transition_cumprobs, transition_targets, transition_step_outs) =
        construct_transitions(
            weights,
            leg_states,
            max_worm_count,
            dims,
            tolerance,
            lp_tolerance,
        )

    return VertexData{NSites}(
        energy_offset,
        dims,
        diagonal_vertices,
        signs,
        weights,
        transitions,
        transition_cumprobs,
        transition_targets,
        transition_step_outs,
        leg_states,
    )
end


vertex_count(vd::VertexData) = length(vd.weights)
get_diagonal_vertex(vd::VertexData, compound_state_idx::Integer) =
    vd.diagonal_vertices[compound_state_idx]
get_sign(vd::VertexData, v::VertexCode) = vd.signs[get_vertex_idx(v)]

function scatter(
    vd::VertexData,
    v::VertexCode,
    leg_in::Integer,
    worm_in::WormIdx,
    random::AbstractFloat,
)
    vi = get_vertex_idx(v)

    t = vd.transitions[leg_in, worm_in, vi]

    @debug begin
        @assert worm_in < worm_count(vd.dims[leg_in%length(vd.dims)])
        @assert !isinvalid(t)
    end

    out =
        findfirst(p -> random < p, @view vd.transition_cumprobs[t.offset:t.offset+t.length])

    (leg_out, worm_out) = vd.transition_step_outs[t.offset+out]

    return (leg_out, worm_out, vd.transition_targets[t.offset+out])
end

get_weight(vd::VertexData, v::VertexCode) = isinvalid(v) ? 0 : vd.weights[get_vertex_idx(v)]


function calc_energy_offset(H::AbstractMatrix, energy_offset_factor::AbstractFloat)
    hmin = minimum(H[collect(diagind(H))])
    hmax = maximum(H[collect(diagind(H))])

    epsilon = (hmax - hmin) * energy_offset_factor
    return -(hmax + epsilon)
end

function construct_vertices(
    dims::NTuple{NSites,<:Integer},
    bond_hamiltonian::AbstractMatrix,
    energy_offset::AbstractFloat,
    tolerance::AbstractFloat,
) where {NSites}
    diagonal_vertices = fill(VertexCode(nothing), size(bond_hamiltonian, 1))
    weights = Vector{eltype(bond_hamiltonian)}()
    leg_states = StateIdx[]
    signs = Int8[]

    for (i, j) in Tuple.(CartesianIndices(bond_hamiltonian))
        w = -bond_hamiltonian[i, j]
        if i == j
            w -= energy_offset
        end

        if abs(w) > tolerance
            if i == j
                diagonal_vertices[i] = VertexCode(true, length(weights))
            end

            for tmp in [i, j]
                tmp -= 1
                for d in dims[end:-1:1]
                    push!(leg_states, tmp % d + 1)
                    tmp ÷= d
                end
                reverse!(@view(leg_states[end-NSites+1:end]))
            end

            push!(weights, abs(w))
            push!(signs, w >= 0 ? 1 : -1)
        end
    end

    return diagonal_vertices, weights, reshape(leg_states, 2 * NSites, :), signs
end

function wrap_vertex_idx(
    leg_states::AbstractMatrix{StateIdx},
    vertex_idx::Union{<:Integer,Nothing},
)
    if vertex_idx !== nothing
        return VertexCode(nothing)
    end
    leg_count = size(leg_states, 1)

    ls = @view leg_states[:, vertex_idx]
    diagonal = true
    for s = 1:leg_count÷2
        diagonal &= ls[s] == ls[leg_count÷2+s]
    end

    return VertexCode(diagonal, vertex_idx)
end

function vertex_apply_change(
    leg_states::AbstractMatrix{StateIdx},
    dims::NTuple{NSites,<:Integer},
    vertex::Integer,
    step_in::Tuple,
    step_out::Tuple,
) where {NSites}
    new_leg_state = leg_states[:, vertex]

    (leg_in, worm_in) = step_in
    (leg_out, worm_out) = step_out

    dim_in = dims[leg_in%NSites]
    dim_out = dims[leg_out%NSites]

    new_leg_state[leg_in] = worm_action(worm_in, new_leg_state[leg_in], dim_in)
    new_leg_state[leg_out] = worm_action(worm_out, new_leg_state[leg_out], dim_out)

    return @views findfirst(v -> new_leg_state == leg_states[:, v], 1:size(leg_states, 2))
end

site_of_leg(leg::Integer, num_sites::Integer) = leg > num_sites ? leg - num_sites : leg

function construct_transitions(
    weights::AbstractArray{<:AbstractFloat},
    leg_states::AbstractMatrix{StateIdx},
    max_worm_count::Integer,
    dims::NTuple{NSites,<:Integer},
    tolerance::AbstractFloat,
    lp_tolerance::AbstractFloat,
) where {NSites}
    used_inaccurate_truncations = false

    leg_count = 2 * NSites

    transitions = Array{Transition,3}(undef, leg_count, max_worm_count, length(weights))
    transition_cumprobs = Float64[]
    transition_targets = VertexCode[]
    transition_step_outs = Tuple{Int,Int}[]

    for array in [transition_cumprobs, transition_targets, transition_step_outs]
        sizehint!(array, length(transitions) * 5) # assume sparseness
    end

    steps = Tuple{Int,Int}[]
    inv_steps = zeros(Int, 2, leg_count, max_worm_count)
    step_idx = -ones(Int, leg_count, max_worm_count)

    for worm = 1:max_worm_count
        for leg = 1:leg_count
            dim = dims[site_of_leg(leg, NSites)]
            if worm < worm_count(dim)
                push!(steps, worm * leg_count + leg)
                step_idx[leg, worm] = length(steps)
                inv_steps[:, leg, worm] .= (leg, worm_inverse(worm, dim))
            end
        end
    end

    variables = NamedTuple{(:step_in, :step_out),Tuple{Int,Int}}[]
    inverse(variable) =
        (step_in = inv_steps[variable.step_out], step_out = inv_steps[variable.step_in])

    for step_in in steps
        for step_out in steps
            vc = (step_in = step_in, step_out = step_out)
            inv_exists = findfirst(x -> x == inverse(x), variables) !== nothing
            if !inv_exists
                push!(variables, vc)
            end
        end
    end

    cost = Float64[vc == inverse(vc) for vc in variables] # discourage bounces

    optimizer = HiGHS.Optimizer()

    x = MOI.add_variables(optimizer, length(variables))
    MOI.set(
        optimizer,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(cost, x), 0.0),
    )
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    #MOI.set(optimizer, MOI.AbsoluteGapTolerance(), lp_tolerance)

    for xi in x
        MOI.add_constraint(optimizer, xi, MOI.GreaterThan(0))
    end

    constraint_idxs = MOI.ConstraintIndex[]
    for step in steps
        xs = Vector{eltype(x)}()
        for (xi, vc) in zip(x, variables)
            if vc.step_in == step || inverse(vc).step_in == step
                push!(xs, xi)
            end
        end

        ci = MOI.add_constraint(
            optimizer,
            MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(ones(length(xs)), xs), 0.0),
            MOI.EqualTo(0.0),
        )
        push!(constraint_idxs, ci)
    end

    targets = Matrix{Int}(undef, leg_count, max_worm_count)
    constraints = zeros(length(steps))

    while true
        v = nothing
        step_in = nothing

        for empty_v in eachindex(weights)
            if v !== nothing
                break
            end

            for empty_step_in in steps
                if isinvalid(transitions[empty_step_in..., empty_v])
                    v = empty_v
                    step_in = empty_step_in
                    break
                end
            end
        end

        if v === nothing # all transitions calculated
            break
        end

        for step_out in steps
            target = vertex_apply_change(leg_states, dims, v, step_in, step_out)
            targets[step_out...] = target
            constraint = target !== nothing ? weights[target] : 0
            MOI.set(
                optimizer,
                MOI.ConstraintSet(),
                constraint_idxs[step_idx[step_out...]],
                MOI.EqualTo(constraint),
            )
            constraints[step_idx[step_out...]] = constraint
        end

        MOI.optimize!(optimizer)

        if !MOI.get(optimizer, MOI.TerminationStatusCode()) == MOI.OPTIMAL
            error(
                "transition probability optimization failed: $(MOI.get(optimizer, MOI.TerminationStatus()))",
            )
        end

        solution = MOI.get(optimizer, MOI.VariablePrimal(), x)

        for in in steps
            if targets[in...] !== nothing
                norm = constraints[step_idx[in...]] == 0 ? 1 : constraints[step_idx[in]]
                offset = length(transition_cumprobs)
                length = 0

                for out in steps
                    var = findfirst(
                        vc ->
                            (vc.step_in == in && vc.step_out == out) || (
                                inverse(vc).step_in == in && inverse(vc).step_out == out
                            ),
                        variables,
                    )
                    @assert var !== nothing

                    prob = solution[var] / norm

                    @assert solution[var] >= -lp_tolerance
                    if prob < 0.0
                        used_inaccurate_truncations = prob < -tolerance
                        prob = 0.0
                    end

                    @assert solution[var] <= norm + lp_tolerance
                    if prob > 1.0
                        used_inaccurate_truncations = prob > 1 + tolerance
                        prob = 1.0
                    end

                    if prob > tolerance / length(steps)
                        push!(transition_cumprobs, prob)

                        in_inv =
                            variables[var].step_in == in ? inverse(variables[var]).step_in :
                            variables[var].step_in
                        push!(
                            transition_targets,
                            wrap_vertex_idx(leg_states, targets[in_inv]),
                        )
                        push!(transition_step_outs, out)
                        @assert !isinvalid(transition_targets[end])

                        length += 1
                    end
                end

                transitions[in..., targets[in...]] = Transition(offset, length)

                probs = @view transition_cumprobs[offset:offset+length]
                cumsum!(probs, probs)
                normed_norm = probs[end]
                if offset >= 0 && abs(normed_norm - 1) > tolerance
                    @error "normalization error: $(normed_norm) ≠ 1"
                end
            end
        end
    end

    if used_inaccurate_truncations
        @warn "had to truncate some probabilities in a possibly inaccurate way!"
    end

    return transitions, transition_cumprobs, transition_targets, transition_step_outs
end
