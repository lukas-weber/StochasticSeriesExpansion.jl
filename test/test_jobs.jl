using Carlo
using Carlo.JobTools

function testjob_magnet_square(sweeps::Integer, thermalization::Integer)
    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 1000

    tm.seed = 124535

    Ts = range(0.04, 4.00, length = 7)

    tm.model = S.MagnetModel
    tm.unitcell = S.UnitCells.square
    tm.Lx = 2
    tm.Ly = 4

    tm.measure = S.all_magnetization_estimators(tm.model, S.dimension(tm.unitcell))

    tm.J = 1.23
    tm.hz = -0.2

    tm.ed_run = 1

    for T in Ts
        task(tm; T = T)
    end

    return make_tasks(tm)
end

function testjob_honeycomb(sweeps::Integer, thermalization::Integer)
    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 1000

    tm.seed = 124535

    Ts = range(0.04, 4.00, length = 7)

    tm.model = S.MagnetModel
    tm.unitcell = S.UnitCells.honeycomb
    tm.Lx = 2
    tm.Ly = 2

    tm.parameter_map = Dict(
        :sites => [Dict(:S => :Sa), Dict(:S => :Sb)],
        :bonds => [Dict(:J => :J1), Dict(:J => :J2), Dict(:J => :J3)],
    )

    tm.measure = S.all_magnetization_estimators(tm.model, S.dimension(tm.unitcell))

    tm.J1 = 1.0
    tm.J2 = 0.5
    tm.J3 = 1.0
    tm.Dz = 0.2
    tm.Dx = 0.5

    tm.Sa = 1 // 2
    tm.Sb = 1

    tm.ed_run = 1

    for T in Ts
        task(tm; T = T)
    end

    make_tasks(tm)
end


function testjob_magnet_bench(sweeps::Integer, thermalization::Integer)
    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 100

    tm.seed = 124535

    tm.T = 0.05

    tm.model = S.MagnetModel
    tm.unitcell = S.UnitCells.square
    tm.Lx = 10
    tm.Ly = 10

    tm.measure = S.all_magnetization_estimators(S.MagnetModel, S.dimension(tm.unitcell))

    tm.J = 1.23
    tm.hz = -0.2

    tm.ed_run = 1

    task(tm)

    return make_tasks(tm)
end

function generate_test_jobs(
    jobdir::AbstractString,
    sweeps::Integer,
    thermalization::Integer,
)
    jobs = Dict{Symbol,JobInfo}()

    function add_test_job(func)
        jobs[Symbol(func)] = JobInfo(
            "$jobdir/$(string(func))",
            S.MC;
            run_time = "15:00",
            checkpoint_time = "15:00",
            tasks = func(sweeps, thermalization),
        )

        return nothing
    end

    add_test_job(testjob_honeycomb)
    add_test_job(testjob_magnet_square)

    return jobs
end
