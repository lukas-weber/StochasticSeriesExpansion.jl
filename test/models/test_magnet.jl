@testset "Magnet" begin
    L = 5
    hz = 0.3
    m = S.MagnetModel(
        Dict(
            :unitcell => S.UnitCells.square,
            :Lx => L,
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
        Dict(:unitcell => S.UnitCells.square, :Lx => L, :Ly => 2 * L, :J => 1.0, :hz => hz),
    )

    @test length(m.site_params) == 2 * L^2
    @test length(m.bond_params) == 4 * L^2
    @test all(m.bond_params[1].hz .≈ (hz / 4, hz / 4))

    Jx = 1.0
    Jy = 2.0
    m = S.MagnetModel(
        Dict(
            :unitcell => S.UnitCells.square,
            :parameter_map =>
                Dict(:bonds => Dict(1 => Dict(:J => :Jx), 2 => Dict(:J => :Jy))),
            :Lx => L,
            :Lx => 2 * L,
            :hz => hz,
            :Jx => Jx,
            :Jy => Jy,
        ),
    )

    @test m.bond_params[1].J ≈ Jx
    @test m.bond_params[2].J ≈ Jy

    @test_throws KeyError S.MagnetModel(
        Dict(
            :unitcell => S.UnitCells.square,
            :parameter_map =>
                Dict(:bonds => Dict(1 => Dict(:J => :Jx), 2 => Dict(:J => :nonexistent))),
            :Lx => L,
            :Lx => 2 * L,
            :hz => hz,
            :Jx => Jx,
            :Jy => Jy,
        ),
    )
end
