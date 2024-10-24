@testset "Magnet" begin
    L = 5
    hz = 0.3
    m = S.MagnetModel(
        Dict(
            :lattice => (; unitcell = S.UnitCells.square, size = (L, L)),
            :J => 1.2,
            :hz => hz,
            :Dz => 0.32,
            :Dx => 0.06,
        ),
    )

    @test length(m.site_params) == L^2
    @test length(m.bond_params) == 2 * L^2
    @test all(m.bond_params[1].hz .≈ (hz / 4, hz / 4))

    m = S.MagnetModel(
        Dict(
            :lattice => (unitcell = S.UnitCells.square, size = (L, 2L)),
            :J => 1.0,
            :hz => hz,
        ),
    )

    @test length(m.site_params) == 2 * L^2
    @test length(m.bond_params) == 4 * L^2
    @test all(m.bond_params[1].hz .≈ (hz / 4, hz / 4))

    Jx = 1.0
    Jy = 2.0
    m = S.MagnetModel(
        Dict(
            :lattice => (; unitcell = S.UnitCells.square, size = (L, 2L)),
            :parameter_map => (J = (:Jx, :Jy),),
            :hz => hz,
            :Jx => Jx,
            :Jy => Jy,
        ),
    )

    @test m.bond_params[1].J ≈ Jx
    @test m.bond_params[2].J ≈ Jy

    @test_throws KeyError S.MagnetModel(
        Dict(
            :lattice => (; unitcell = S.UnitCells.square, size = (L, 2L)),
            :parameter_map => (; J = (:Jx, :nonexistent),),
            :hz => hz,
            :Jx => Jx,
            :Jy => Jy,
        ),
    )
end
