using LoadLeveller.ResultTools
using Distributions
using HypothesisTests
using Measurements

include("ed/ed.jl")
include("ed/magnet.jl")
include("test_jobs.jl")


function run_mc(job::JobInfo)
    LoadLeveller.start(LoadLeveller.SingleRunner{job.mc}, job)

    return DataFrame(ResultTools.dataframe("$(job.dir)/../$(job.name).results.json"))
end

@testset "ED Comparison" begin
    mktempdir() do tmpdir
        jobs = generate_test_jobs(tmpdir, 100000, 10000)

        for (jobname, job) in jobs
            @testset "$(job.name)" begin
                get_model(::Type{S.MC{Model,NSites}}) where {Model,NSites} = Model
                model = get_model(job.mc)
                (obsnames, ed_data) = run_ed(model, job)
                mc_data = run_mc(job)

                for obsname in obsnames
                    # hack for the case that the error is exactly zero
                    mc_data_nudge =
                        [m.val ± (m.err != 0 ? m.err : 1e-9) for m in mc_data[!, obsname]]

                    z = stdscore.(mc_data_nudge, ed_data[!, obsname])
                    p = pvalue(ExactOneSampleKSTest(z, Normal()))
                    if p <= 1e-3
                        println("$obsname: p = $p")
                        println("MC: $(mc_data_nudge)\nED:$(ed_data[!,obsname])")
                        gui(plot(ed_data[!,:T], [mc_data_nudge, ed_data[!,obsname]]))
                    end

                    @test p > 1e-3
                end
            end
        end
    end
end
