@testset "ClusterModel" begin

    @testset "lift_twobody_operator" begin
        A = reshape(1:16, 4, 4)
        B = reshape(1:9, 3, 3)
        C = reshape(1:4, 2, 2)

        dims = size.((A, B, C), 1)

        @test S.lift_twobody_operator(kron(A, B), dims, (1, 2)) == kron(A, B, I(dims[3]))
        @test S.lift_twobody_operator(kron(A, B), size.((B, C, A), 1), (3, 1)) *
              kron(I(dims[2]), C, I(dims[1])) == kron(B, C, A)
        @test S.lift_twobody_operator(kron(A, B), size.((B, C, A), 1), (3, 1)) *
              kron(B, C, I(dims[1])) == kron(B^2, C, A)
    end

    @testset "dimer magnet" begin
        params = Dict(
            :inner_model => S.MagnetModel,
            :lattice =>
                (unitcell = S.UnitCells.fully_frust_square_bilayer, size = (2, 2)),
            :parameter_map => (;
                J = [:J1, :J2, :J3, :J4, :J5, :J6, :J7, :J8, :J9],
                Dz = [:Dz1, :Dz2],
            ),
            :cluster_bases => (S.ClusterBases.dimer,),
            :measure_quantum_numbers => [],
            :J1 => 1,
            :J2 => 1.1,
            :J3 => 1.2,
            :J4 => 1.3,
            :J5 => 1.4,
            :J6 => 1.5,
            :J7 => 1.6,
            :J8 => 1.7,
            :J9 => 1.8,
        )

        model = S.ClusterModel(params)

        intracluster_hamiltonians, intercluster_hamiltonians = S.build_cluster_hamiltonians(
            model.cluster_ids[eachindex(model.inner_model.lattice.uc.sites)],
            model.inner_model.lattice.uc.bonds,
            model.inner_model.bond_params,
            model.inner_model.site_params,
        )

        U = model.basis[1].transformation
        @test U' * intracluster_hamiltonians[1] * U ≈
              Diagonal([-3 // 4, 1 // 4, 1 // 4, 1 // 4])

        @test Set(keys(intercluster_hamiltonians)) ==
              Set([(iuc = 1, juc = 1, jd = (1, 0)), (iuc = 1, juc = 1, jd = (0, 1))])

        Hinter = map(h -> kron(U, U)' * h * kron(U, U), values(intercluster_hamiltonians))

        for H in Hinter
            H = reshape(H, 4, 4, 4, 4)
            ε = 1e-10
            H[abs.(H).<ε] .= 0

            Htriplet = reshape(H[2:4, 2:4, 2:4, 2:4], 9, 9)
            Htriplet ./= Htriplet[1, 1]
            @test Htriplet ≈ [
                1 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 0 0 0
                0 0 -1 0 1 0 0 0 0
                0 1 0 0 0 0 0 0 0
                0 0 1 0 0 0 1 0 0
                0 0 0 0 0 0 0 1 0
                0 0 0 0 1 0 -1 0 0
                0 0 0 0 0 1 0 0 0
                0 0 0 0 0 0 0 0 1
            ]

            @test iszero(H[1, :, 1, :])
            @test isdiag(H[:, 1, 1, :])
            @test iszero(H[:, 1, :, 1])
            @test isdiag(H[1, :, :, 1])
        end

        sse_data = S.generate_sse_data(model)
        @test length(sse_data.bonds) == 8
    end
end
