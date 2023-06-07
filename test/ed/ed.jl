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
    for i = 1:size(en.psi, 2)
        res += (en.psi[:, i]' * A * en.psi[:, i]) * en.ρ[i, :]
    end

    if all(imag.(res) .≈ 0.0)
        return real.(res)
    end
    return res
end

function diag_mean(en::Ensemble, Adiag::AbstractVector)
    return vec(sum(Adiag .* en.ρ, dims = 1))
end

"""
    integrated_correlator(ens, A, B)

Compute the integrated correlation function

```math
\\int^\\beta_0 d\\tau \\langle A(\\tau) B\\rangle.
```
"""
function integrated_correlator(
    ens::Ensemble{F},
    A::AbstractMatrix,
    B::AbstractMatrix,
) where {F}
    Anm = ens.psi' * A * ens.psi
    Bnm = ens.psi' * B * ens.psi

    return [
        sum(
            2ens.ρ[n, i] / (ens.Es[n] == ens.Es[m] ? 2T : ens.Es[m] - ens.Es[n]) *
            Anm[n, m] *
            Bnm[m, n] for n in eachindex(ens.Es), m in eachindex(ens.Es)
        ) for (i, T) in enumerate(ens.Ts)
    ]
end

@struct_equal S.Models.Magnet

function summarize_tasks(job::JobInfo)
    summarized_tasks = Tuple{Vector{String},Vector{Float64},Dict{Symbol,Any}}[]

    for task in job.tasks

        params = deepcopy(task.params)
        T = pop!(params, :T)

        existing = findfirst(t -> t[3][:ed_run] == params[:ed_run], summarized_tasks)
        if existing === nothing
            push!(summarized_tasks, ([task.name], [T], params))
        else
            push!(summarized_tasks[existing][1], task.name)
            push!(summarized_tasks[existing][2], T)
        end
    end

    return summarized_tasks
end


function run_ed(::Type{Model}, job::JobInfo) where {Model}
    results = Dict{Symbol,Any}[]

    obsnames = Set{Symbol}()

    for (names, Ts, params) in summarize_tasks(job)
        model = Model(params)
        H = hamiltonian(model)

        (Es, psi) = eigen(Symmetric(Matrix(H)))

        ensemble = Ensemble(Ts, Es, psi)

        energy = diag_mean(ensemble, Es) / S.normalization_site_count(model)
        specific_heat =
            (diag_mean(ensemble, Es .^ 2) - diag_mean(ensemble, Es) .^ 2) ./
            (Ts .^ 2 * S.normalization_site_count(model))

        obs = Dict(:Energy => energy, :SpecificHeat => specific_heat)

        calc_observables!(obs, model, ensemble)

        union!(obsnames, keys(obs))
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

    return collect(obsnames), DataFrame(results)
end
