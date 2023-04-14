using LoadLeveller.JobTools
using LinearAlgebra
using SparseArrays
using StructEquality
using DataFrames

struct Lifter
    site_dims::Vector{Int}
    cum_site_dims::Vector{Int}
end

Lifter(site_dims::AbstractVector{<:Integer}) =
    Lifter(site_dims, vcat([1], cumprod(site_dims)))

function lift(l::Lifter, pos::Integer, op::AbstractMatrix)
    @assert all(size(op) .== l.site_dims[pos])
    dleft = l.cum_site_dims[pos]
    dright = l.cum_site_dims[end] ÷ l.cum_site_dims[pos+1]
    return kron(kron(sparse(I, dleft, dleft), op), sparse(I, dright, dright))
end

id(l::Lifter) = sparse(I, l.cum_site_dims[end], l.cum_site_dims[end])

function spin(l::Lifter, pos::Integer, idx::Integer)
    (splus, sz) = S.spin_operators(l.site_dims[pos])

    if idx == 1
        return lift(l, pos, 0.5 * sparse(splus + splus'))
    elseif idx == 2
        return lift(l, pos, 0.5im * sparse(splus - splus'))
    elseif idx == 3
        return lift(l, pos, sparse(sz))
    end
    error("invalid spin index")
end

function heisen_bond(l::Lifter, posi::Integer, posj::Integer)
    return real.(sum([spin(l, posi, a) * spin(l, posj, a) for a = 1:3]))
end


function hamiltonian(magnet::S.Models.Magnet)
    dims = [Int(site.spin_mag * 2 + 1) for site in magnet.site_params]

    lifter = Lifter(dims)
    H = spzeros(prod(dims), prod(dims))
    for (bond, params) in zip(magnet.lattice.bonds, magnet.bond_params)
        H +=
            params.J * heisen_bond(lifter, bond.i, bond.j) +
            params.d * spin(lifter, bond.i, 3) * spin(lifter, bond.j, 3) +
            params.hz * (spin(lifter, bond.i, 3) + spin(lifter, bond.j, 3))
        +params.Dx * (spin(lifter, bond.i, 1)^2 + spin(lifter, bond.j, 1)^2) +
        params.Dz * (spin(lifter, bond.i, 1)^2 + spin(lifter, bond.j, 3)^2)
    end

    return H
end


struct Ensemble{F,CF}
    Ts::Vector{F}
    Es::Vector{F}
    psi::Matrix{CF}
    ρ::Matrix{F}
end

function Ensemble(Ts::AbstractVector, Es::AbstractVector, psi::AbstractMatrix)
    Enorm = Es .- minimum(Es)
    ρ = exp.(-Enorm ./ Ts')
    Z = sum(ρ, dims = 1)
    ρ ./= Z

    return Ensemble(Ts, Es, psi, ρ)
end

function mean(en::Ensemble{F}, A::AbstractMatrix) where {F}
    res = zeros(Complex{F}, length(en.Ts))
    for i in size(en.psi, 2)
        res += en.psi[:, i]' * A * en.psi[:, i] * en.ρ[i, :]
    end

    if all(imag.(res) ≈ 0)
        return real.(res)
    end
    return res
end

function diag_mean(en::Ensemble, Adiag::AbstractVector)
    return sum(Adiag .* en.ρ, dims = 1)
end

@struct_equal S.Models.Magnet

function summarize_tasks(job::JobInfo)
    summarized_tasks = Tuple{Vector{String},Vector{Float64},Dict{Symbol,Any}}[]

    for task in job.tasks

        params = deepcopy(task.params)
        T = pop!(params, :T)

        existing = findfirst(t -> t[2] == params, summarized_tasks)
        if existing === nothing
            push!(summarized_tasks, ([task.name], [T], params))
        else
            push!(summarized_task[1], task.name)
            push!(summarized_task[2], T)
        end
    end

    return summarized_tasks
end


function run_ed(::Type{Model}, job::JobInfo) where {Model}
    results = Dict{Symbol,Any}[]

    for (names, Ts, params) in summarize_tasks(job)
        model = Model(params)
        H = hamiltonian(model)

        (Es, psi) = eigen(Symmetric(Matrix(H)))

        ensemble = Ensemble(Ts, Es, psi)

        Emean = diag_mean(ensemble, Es) / S.normalization_site_count(model)

        obs = Dict(:Energy => Emean)

        for (i, T) in enumerate(Ts)
            push!(
                results,
                merge(
                    Dict(:task => names[i], :T => T),
                    params,
                    Dict(obsname => obsval[i] for (obsname, obsval) in obs),
                ),
            )
        end
    end

    return DataFrame(results)
end
