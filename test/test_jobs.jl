using LoadLeveller
using LoadLeveller.JobTools
using LoadLeveller.ResultTools
using Distributions
using HypothesisTests
using Measurements

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
            S.mc(Model);
            run_time = "15:00",
            checkpoint_time = "15:00",
            tasks = make_tasks(tm),
        )
    end

    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 1000

    Ts = range(0.02, 4.00, length = 7)

    tm.unitcell = S.UnitCells.square
    tm.Lx = 2
    tm.Ly = 4

    tm.J = 1.23
    tm.hz = -0.4

    for T in Ts
        task(tm; T = T)
    end

    push!(jobs, test_job("magnet_square", S.Models.Magnet, tm))

    return jobs
end

function run_mc(job::JobInfo)
    LoadLeveller.start(LoadLeveller.SingleRunner{job.mc}, job)

    return DataFrame(ResultTools.dataframe("$(job.dir)/../$(job.name).results.json"))
end

@testset "ED Comparison" begin
    mktempdir() do tmpdir
        jobs = generate_test_jobs(tmpdir, 100000, 10000)

        for job in jobs
            @testset "$(job.name)" begin
                get_model(::Type{S.MC{Model,NSites}}) where {Model,NSites} = Model
                model = get_model(job.mc)
                (obsnames, ed_data) = run_ed(model, job)
                mc_data = run_mc(job)

                for obsname in obsnames
                    z = stdscore.(mc_data[!, obsname], ed_data[!, obsname])
                    p = pvalue(ExactOneSampleKSTest(z, Normal()))
                    if p <= 1e-3
                        println("$obsname: p = $p")
                    end

                    @test p > 1e-3
                end
            end
        end
    end
end
