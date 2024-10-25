using Carlo.ResultTools
using Measurements

include("ed/ed.jl")
include("ed/magnet.jl")


function run_mc(job::JobInfo)
    Carlo.start(Carlo.SingleScheduler, job)

    return DataFrame(ResultTools.dataframe("$(job.dir)/../$(job.name).results.json"))
end

@testset "ED Comparison" begin
    mktempdir() do tmpdir
        jobs = generate_test_jobs(tmpdir, 40000, 5000)

        for (jobname, job) in jobs
            @testset "$(job.name)" begin
                println("Running ED...")
                (obsnames, ed_data) = run_ed(job)
                println("Running MC...")
                mc_data = run_mc(job)

                for obsname in obsnames
                    # hack for the case that the error is exactly zero
                    mc_data_nudge =
                        [m.val ± (m.err != 0 ? m.err : 1e-8) for m in mc_data[!, obsname]]

                    z = abs.(stdscore.(mc_data_nudge, ed_data[!, obsname]))
                    if any(>(4), z)
                        println("$obsname:")
                        display(
                            DataFrame(
                                :MC => mc_data_nudge,
                                :ED => ed_data[!, obsname],
                                :z => z,
                            ),
                        )
                        if isdefined(Main, :Plots)
                            pl = plot(
                                ed_data[!, :T],
                                ed_data[!, obsname],
                                label = "$obsname: ED",
                            )
                            plot!(
                                pl,
                                ed_data[!, :T],
                                mc_data_nudge,
                                label = "$obsname: QMC",
                            )
                            gui(pl)
                        end
                    end

                    @test all(<=(4), z)
                end
            end
        end
    end
end
