
@testset "SSEData" begin
    @testset "inconsistent dims" begin
        vd32 = S.VertexData((3, 2), Matrix(1.0 * I, 3 * 2, 3 * 2))
        vd23 = S.VertexData((2, 3), Matrix(1.0 * I, 3 * 2, 3 * 2))

        vds = [vd32, vd23]

        bonds_correct = [S.SSEBond(1, (1, 2)), S.SSEBond(2, (2, 3))]
        bonds_incorrect1 = [S.SSEBond(1, (1, 2)), S.SSEBond(1, (2, 3))]
        bonds_incorrect2 = [S.SSEBond(1, (1, 2)), S.SSEBond(2, (3, 2))]

        sites = S.generate_sites_from_bonds(vds, bonds_correct)

        @test length(sites) == 3
        @test [s.dim for s in sites] == [3, 2, 3]

        @test_throws ErrorException S.generate_sites_from_bonds(vds, bonds_incorrect1)
        @test_throws ErrorException S.generate_sites_from_bonds(vds, bonds_incorrect2)
    end
end
