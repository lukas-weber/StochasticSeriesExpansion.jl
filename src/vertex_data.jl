using HiGHS
using LinearAlgebra
using Printf
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

"""
    VertexData(dims::Tuple{Integer...}, bond_hamiltonian::AbstractMatrix; 
               energy_offset_factor = 0.25, tolerance = 1e-7, lp_tolerance = 1e-10)

This object holds the probability tables used for a bond in the abstract loop update algorithm.

`dims` are the local Hilbert space dimensions of each site of the bond.

`bond_hamiltonian` is the `prod(dims)`×`prod(dims)` matrix describing the bond Hamiltonian in the computational basis. The representation of the product Hilbert space follows the convention of `LinearAlgebra.kron`, so that if you have site-local operators `op1`, `op2`, `kron(op1, op2)` will give you a correct bond Hamiltonian.

A constant is added to the bond hamiltonian so that its smallest diagonal element is
```math
\\min_i(h_{ii}) = \\alpha \\times [\\max_i (h_{ii}) - \\min_i (h_{ii})] \\ge 0,
```
where ``\\alpha`` is given by the parameter `energy_offset_factor`.

Matrix elements smaller `tolerance` are considered zero. `lp_tolerance` sets the tolerance for the underlying linear programming problem.
"""
function VertexData(
    dims::NTuple{NSites,<:Integer},
    bond_hamiltonian::AbstractMatrix;
    energy_offset_factor::AbstractFloat = 0.25,
    tolerance::AbstractFloat = 1e-7,
    lp_tolerance::AbstractFloat = 1e-10,
) where {NSites}
    energy_offset = calc_energy_offset(bond_hamiltonian, energy_offset_factor)
    total_dim = prod(dims)
    max_worm_count = maximum(worm_count, dims)

    (diagonal_vertices, weights, leg_states, signs) =
        construct_vertices(dims, bond_hamiltonian, energy_offset, tolerance)

    @assert NSites >= 1
    @assert total_dim == size(bond_hamiltonian, 1) == size(bond_hamiltonian, 2)
    @assert size(leg_states) == (2 * NSites, length(weights))

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
get_vertex_weight(vd::VertexData, v::VertexCode) =
    isinvalid(v) ? 0.0 : vd.weights[get_vertex_idx(v)]
get_sign(vd::VertexData, v::VertexCode) = vd.signs[get_vertex_idx(v)]

"""
    get_leg_state(vd::VertexData, v::VertexCode)

Returns the leg states (matrix elements) of the vertex `v` as an array of state indices.
"""
get_leg_state(vd::VertexData, v::VertexCode) = @view vd.leg_states[:, get_vertex_idx(v)]

function scatter(
    vd::VertexData{NSites},
    v::VertexCode,
    leg_in::Integer,
    worm_in::WormIdx,
    random::AbstractFloat,
) where {NSites}
    vi = get_vertex_idx(v)

    @inbounds t = vd.transitions[leg_in, worm_in, vi]

    @inbounds for out = t.offset:t.offset+t.length
        if random < vd.transition_cumprobs[out]
            (leg_out, worm_out) = vd.transition_step_outs[out]

            return (leg_out, worm_out, vd.transition_targets[out])
        end
    end
    return (-1, -1, VertexCode(nothing))
end



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
                reversed_i = join_idx(dims, reverse(split_idx(reverse(dims), i)))
                diagonal_vertices[reversed_i] = VertexCode(true, length(weights) + 1)
            end

            append!(leg_states, reverse(split_idx(reverse(dims), i)))
            append!(leg_states, reverse(split_idx(reverse(dims), j)))

            push!(weights, abs(w))
            push!(signs, w >= 0 ? 1 : -1)
        end
    end

    return diagonal_vertices, weights, reshape(leg_states, 2NSites, :), signs
end

function wrap_vertex_idx(
    leg_states::AbstractMatrix{StateIdx},
    vertex_idx::Union{<:Integer,Nothing},
)
    if vertex_idx === nothing
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


site_of_leg(leg::Integer, num_sites::Integer) = leg > num_sites ? leg - num_sites : leg

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

    dim_in = dims[site_of_leg(leg_in, NSites)]
    dim_out = dims[site_of_leg(leg_out, NSites)]

    new_leg_state[leg_in] = worm_action(worm_in, new_leg_state[leg_in], dim_in)
    new_leg_state[leg_out] = worm_action(worm_out, new_leg_state[leg_out], dim_out)

    return @views findfirst(v -> new_leg_state == leg_states[:, v], 1:size(leg_states, 2))
end

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

    transitions = fill(Transition(), leg_count, max_worm_count, length(weights))
    transition_cumprobs = Float64[]
    transition_targets = VertexCode[]
    transition_step_outs = Tuple{Int,Int}[]

    for array in [transition_cumprobs, transition_targets, transition_step_outs]
        sizehint!(array, length(transitions) * 5) # assume sparseness
    end

    steps = Tuple{Int,Int}[]
    inv_steps = fill((0, 0), leg_count, max_worm_count)
    step_idx = -ones(Int, leg_count, max_worm_count)

    for worm = 1:max_worm_count
        for leg = 1:leg_count
            dim = dims[site_of_leg(leg, NSites)]
            if worm <= worm_count(dim)
                push!(steps, (leg, worm))
                step_idx[leg, worm] = length(steps)
                inv_steps[leg, worm] = (leg, worm_inverse(worm, dim))
            end
        end
    end

    variables = NamedTuple{(:step_in, :step_out),Tuple{Tuple{Int,Int},Tuple{Int,Int}}}[]
    inverse(variable) = (
        step_in = inv_steps[variable.step_out...],
        step_out = inv_steps[variable.step_in...],
    )

    for step_in in steps
        for step_out in steps
            vc = (step_in = step_in, step_out = step_out)
            inv_exists = findfirst(x -> x == inverse(vc), variables) !== nothing
            if !inv_exists
                push!(variables, vc)
            end
        end
    end

    cost = Float64[vc == inverse(vc) for vc in variables] # discourage bounces

    optimizer = HiGHS.Optimizer()

    MOI.set(optimizer, MOI.Silent(), true)
    x = MOI.add_variables(optimizer, length(variables))
    MOI.set(
        optimizer,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(cost, x), 0.0),
    )
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("primal_feasibility_tolerance"),
        lp_tolerance,
    )
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("dual_feasibility_tolerance"),
        lp_tolerance,
    )

    for xi in x
        MOI.add_constraint(optimizer, xi, MOI.GreaterThan(0.0))
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

    targets = Matrix{Union{Nothing,Int}}(undef, leg_count, max_worm_count)
    constraints = zeros(length(steps))

    while true
        v = nothing
        step_in = nothing

        for empty_v = 1:length(weights)
            iempty_step_in = findfirst(s -> isinvalid(transitions[s..., empty_v]), steps)
            if iempty_step_in !== nothing
                v = empty_v
                step_in = inv_steps[steps[iempty_step_in]...]
                break
            end
        end

        if v === nothing
            break # all transitions calculated
        end

        for step_out in steps
            target = vertex_apply_change(leg_states, dims, v, step_in, step_out)
            targets[step_out...] = target
            constraint = target !== nothing ? weights[target] : 0.0
            MOI.set(
                optimizer,
                MOI.ConstraintSet(),
                constraint_idxs[step_idx[step_out...]],
                MOI.EqualTo(constraint),
            )
            constraints[step_idx[step_out...]] = constraint
        end

        MOI.optimize!(optimizer)

        if MOI.get(optimizer, MOI.TerminationStatus()) != MOI.OPTIMAL
            error(
                "transition probability optimization failed: $(MOI.get(optimizer, MOI.TerminationStatus()))",
            )
        end

        solution = MOI.get(optimizer, MOI.VariablePrimal(), x)

        for in in steps
            if targets[in...] !== nothing
                norm = constraints[step_idx[in...]] == 0 ? 1 : constraints[step_idx[in...]]
                offset = length(transition_cumprobs) + 1
                len = -1 # one-based indexing...

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
                            wrap_vertex_idx(leg_states, targets[in_inv...]),
                        )
                        push!(transition_step_outs, out)
                        @assert !isinvalid(transition_targets[end])

                        len += 1
                    end
                end

                transitions[in..., targets[in...]] = Transition(offset, len)

                probs = @view transition_cumprobs[offset:offset+len]
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

function Base.show(io::IO, vd::VertexData{NSites}) where {NSites}
    println(io, "VertexData(")
    println(io, "weights = $(vd.weights)")
    for v = 1:size(vd.transitions, 3)
        println("vertex $(v): [$(join(string.(vd.leg_states[:,v]), ","))]")
        for worm = 1:size(vd.transitions, 2)
            for leg_in = 1:size(vd.transitions, 1)
                t = vd.transitions[leg_in, worm, v]
                if isinvalid(t)
                    continue
                end

                probs = diff(vd.transition_cumprobs[t.offset:t.offset+t.length])
                targets = vd.transition_targets[t.offset:t.offset+t.length]
                step_outs = vd.transition_step_outs[t.offset:t.offset+t.length]

                @printf(
                    io,
                    "(%d,%d) %s | %s\n",
                    leg_in,
                    worm,
                    join([@sprintf("%.2f", p) for p in probs], " "),
                    join(
                        [
                            @sprintf("(%d,%d,%d)", sout..., get_vertex_idx(tar)) for
                            (sout, tar) in zip(step_outs, targets)
                        ],
                        " ",
                    ),
                )
            end
        end
    end
    print(io, ")")
end
