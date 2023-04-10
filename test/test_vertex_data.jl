
@testset "VertexData" begin
    @testset "S=1/2 Heisenberg" begin
        (splus, sz) = S.spin_operators(2)
        Hbond = kron(sz, sz) + 0.5 * (kron(splus, splus') + kron(splus', splus))

        vd = S.VertexData((2, 2), Hbond; energy_offset_factor = 0.0)

        @show vd.diagonal_vertices
        @test vd.energy_offset ≈ -0.25
        @test sum(S.isinvalid.(vd.diagonal_vertices)) == 2
        @test all(vd.weights .≈ 0.5)

        @show vd.transition_cumprobs
    end

end
