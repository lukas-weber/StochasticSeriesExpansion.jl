import StochasticSeriesExpansion as S

include("test_jobs.jl")

function bench_vertex_list()
    job = generate_test_jobs("", 10000, 10000)["magnet_square"]

    mag = S.MagnetModel(job.tasks[1].params)

    sse_data = S.generate_sse_data(mag)

    v = S.VertexCode(false, 1)
    operators = [S.OperCode(rand(1:length(sse_data.bonds)), v) for _ = 1:1000000]

    vl = S.VertexList{2}(length(sse_data.sites))
    return () -> S.make_vertex_list!(vl, operators, sse_data.bonds)
end

function bench_total()
    function bench()
        mktempdir() do dir
            name, model, tm = testjob_magnet_bench(1000, 1000)
            job = JobInfo(
                "$dir/$name",
                S.MC,
                tasks = make_tasks(tm),
                checkpoint_time = "15:00",
                run_time = "15:00",
            )

            Carlo.start(Carlo.SingleScheduler, job)
        end
        return nothing
    end

    return bench
end
