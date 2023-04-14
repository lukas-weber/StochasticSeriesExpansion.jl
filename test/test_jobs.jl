using LoadLeveller
using LoadLeveller.JobTools
using LoadLeveller.ResultTools

include("ed/ed.jl")

function generate_test_jobs(
    jobdir::AbstractString,
    sweeps::Integer,
    thermalization::Integer,
)
    jobs = JobInfo[]

    function test_job(name::AbstractString, ::Type{Model}, tm::TaskMaker) where {Model}
        return JobInfo(
            "$jobdir/$name",
            S.MC{Model};
            run_time = "15:00",
            checkpoint_time = "15:00",
            tasks = make_tasks(tm),
        )
    end

    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 100

    Ts = range(0.02, 1.0, length = 7)

    tm.unitcell = S.UnitCells.square
    tm.Lx = 2
    tm.J = 1.24
    tm.hz = 0.1

    for T in Ts
        task(tm; T = T)
    end

    push!(jobs, test_job("magnet_square", S.Models.Magnet, tm))

    return jobs
end

function run_mc(job::JobInfo)
    LoadLeveller.start(LoadLeveller.SingleRunner{job.mc}, job)

    return ResultTools.dataframe(job.jobdir * "/results.json")
end

@testset "ED Comparison" begin
    mktempdir() do tmpdir
        jobs = generate_test_jobs(tmpdir, 10000, 4000)

        for job in jobs
            @testset "$(job.name)" begin
                get_model(::Type{S.MC{Model}}) where {Model} = Model
                model = get_model(job.mc)
                ed_data = run_ed(model, job)
                mc_data = run_mc(job)
                @show ed_data, mc_data
            end
        end
    end
end
