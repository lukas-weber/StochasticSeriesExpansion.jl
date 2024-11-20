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
    tm.lattice = (unitcell = S.UnitCells.square, size = (2, 4))

    tm.measure = S.all_magnetization_estimators(tm.model, S.dimension(tm.lattice.unitcell))

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
    tm.lattice = (unitcell = S.UnitCells.honeycomb, size = (2, 2))

    tm.parameter_map = (S = [:Sa, :Sb], J = [:J1, :J2, :J3])

    tm.measure = S.all_magnetization_estimators(tm.model, S.dimension(tm.lattice.unitcell))

    tm.J1 = 1.0
    tm.J2 = 0.5
    tm.J3 = 1.0

    tm.d = 0.2

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

function testjob_fully_frustrated_bilayer(sweeps::Integer, thermalization::Integer)
    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 1000

    tm.seed = 124535

    Ts = range(0.04, 4.00, length = 7)

    tm.model = S.ClusterModel
    tm.inner_model = S.MagnetModel
    tm.cluster_bases = (S.ClusterBases.dimer,)
    tm.lattice = (unitcell = S.UnitCells.fully_frust_square_bilayer, size = (2, 2))

    tm.measure_quantum_numbers = [(name = Symbol(), quantum_number = 2)]

    tm.parameter_map =
        (S = [:Sa, :Sb], J = [:JD, :JxP, :JxP, :JyP, :JyP, :JxX, :JxX, :JyX, :JyX])

    tm.JD = 1
    tm.JxP = 0.55
    tm.JyP = 0.6
    tm.JxX = 0.7
    tm.JyX = 0.3

    tm.d = -0.2

    tm.Dz = 0.2
    tm.Dx = 0.5

    tm.Sa = 1 // 2
    tm.Sb = 1 // 2

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
    tm.lattice = (unitcell = S.UnitCells.square, size = (10, 10))

    tm.measure =
        S.all_magnetization_estimators(S.MagnetModel, S.dimension(tm.lattice.unitcell))

    tm.J = 1.23
    tm.hz = -0.2

    tm.ed_run = 1

    task(tm)

    return make_tasks(tm)
end

function testjob_cluster_bench(sweeps::Integer, thermalization::Integer)
    tm = TaskMaker()
    tm.sweeps = sweeps
    tm.thermalization = thermalization
    tm.binsize = 100

    tm.seed = 124535

    tm.T = 0.1

    tm.model = S.ClusterModel
    tm.inner_model = S.MagnetModel
    tm.cluster_bases = (S.ClusterBases.dimer,)
    tm.lattice = (unitcell = S.UnitCells.fully_frust_square_bilayer, size = (30, 30))

    tm.measure_quantum_numbers = [(name = Symbol(), quantum_number = 2)]

    tm.parameter_map = (S = [:Sa, :Sb], J = vcat([:JD], repeat([:JP], 8)))

    tm.JD = 0.5
    tm.JP = 1

    tm.Sa = 1 // 2
    tm.Sb = 1 // 2

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
    add_test_job(testjob_fully_frustrated_bilayer)

    return jobs
end
