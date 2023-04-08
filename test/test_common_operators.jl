using LinearAlgebra

@testset "Common operators" begin
    (σp, σz) = S.spin_operators(2)
    @test σp ≈ [0 1; 0 0]
    @test σz ≈ [0.5 0; 0 -0.5]

    dims = rand(2:100, 10)

    for d in dims
        (splus, sz) = S.spin_operators(d)
        s = (d - 1) / 2

        @test sz * sz + 0.5 * (splus' * splus + splus * splus') ≈
              diagm(fill(s * (s + 1), d))

        a = S.bosonic_a_operator(d)

        @test a' * a ≈ diagm(0:d-1)
    end
end
