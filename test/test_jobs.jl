using LoadLeveller
using LoadLeveller.JobTools

function testjob_magnet_square(sweeps::Integer, thermalization::Integer)
    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 1000

    Ts = range(0.04, 4.00, length = 7)

    tm.unitcell = S.UnitCells.square
    tm.Lx = 2
    tm.Ly = 4

    tm.measure = S.all_magests(S.Models.Magnet, S.dimension(tm.unitcell))

    tm.J = 1.23
    tm.hz = -0.2

    tm.ed_run = 1

    for T in Ts
        task(tm; T = T)
    end

    return tm
end

function generate_test_jobs(
    jobdir::AbstractString,
    sweeps::Integer,
    thermalization::Integer,
)
    jobs = Dict{String,JobInfo}()

    function add_test_job(
        name::AbstractString,
        model::Type{<:S.AbstractModel},
        tm::TaskMaker,
    )
        jobs[name] = JobInfo(
            "$jobdir/$name",
            S.mc(model);
            run_time = "15:00",
            checkpoint_time = "15:00",
            tasks = make_tasks(tm),
        )

        return nothing
    end

    add_test_job(
        "magnet_square",
        S.Models.Magnet,
        testjob_magnet_square(sweeps, thermalization),
    )

    return jobs
end
