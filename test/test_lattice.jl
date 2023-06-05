using LinearAlgebra

@testset "UnitCell" begin
    uc = S.UnitCells.square
    @test uc.sites[1].sublattice_sign == 1
    @test uc.sites[1].coordination == 4

    uc = S.UnitCells.columnar_dimer
    @test uc.sites[1].sublattice_sign == 1
    @test uc.sites[2].sublattice_sign == -1
    @test uc.sites[1].coordination == 4
    @test uc.sites[2].coordination == 4

    uc = S.UnitCells.honeycomb
    @test uc.sites[1].sublattice_sign == 1
    @test uc.sites[2].sublattice_sign == -1
    @test uc.sites[1].coordination == 3
    @test uc.sites[2].coordination == 3
end

function adjacency_matrix(l::S.Lattice)
    A = zeros(Int, S.site_count(l), S.site_count(l))

    for b in l.bonds
        A[b.i, b.j] = 1
        A[b.j, b.i] = 1
    end

    return A
end

# thanks to https://planetmath.org/enumerationoflatticewalks
num_square_walks(n::Integer) = iseven(n) ? binomial(n, n ÷ 2)^2 : 0
num_honeycomb_walks(n::Integer) =
    iseven(n) ? sum(binomial(n ÷ 2, k)^2 * binomial(2k, k) for k = 0:n÷2) : 0
num_triangle_walks(n::Integer) =
    sum((-2)^(n - k) * binomial(n, k) * sum(binomial(k, j)^3 for j = 0:k) for k = 0:n)



@testset "Lattice" begin

    @testset "n-walks" begin
        for Ls in [(3, 4), (8, 4), (5, 5), (8, 8)]
            ucs = [
                S.UnitCells.square,
                S.UnitCells.columnar_dimer,
                S.UnitCells.honeycomb,
                S.UnitCells.triangle,
            ]
            num_walks = [
                num_square_walks,
                num_square_walks,
                num_honeycomb_walks,
                num_triangle_walks,
            ]

            lattices = (uc -> S.Lattice(uc, Ls)).(ucs)
            As = adjacency_matrix.(lattices)
            (Asquare, Acolumnar_dimer, Ahoneycomb, Atriangle) = As

            for (A, lattice) in zip(As, lattices)
                @test tr(A^2) / 2 == length(lattice.bonds)
            end


            for n = 2:6
                if any(Ls .<= n)
                    continue
                end

                for (lattice, A, nwalks) in zip(lattices, As, num_walks)
                    @test tr(A^n) / S.site_count(lattice) == nwalks(n)
                end
            end
        end
    end

    @testset "neel_vector" begin
        @test S.neel_vector(S.UnitCells.square) == ((1, 1), 0)
        @test S.neel_vector(S.UnitCells.columnar_dimer) == ((1, 0), 1)
        @test S.neel_vector(S.UnitCells.honeycomb) == ((0, 0), 1)
        @test S.neel_vector(S.UnitCells.triangle) === nothing
    end
end
