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

@testset "Lattice" begin

    @testset "n-walks" begin
        for Ls in [(3,4), (8,4), (5,5), (8,8)]
            ucs = [S.UnitCells.square, S.UnitCells.columnar_dimer, S.UnitCells.honeycomb]
            lattices = (uc->S.Lattice(uc, Ls)).(ucs)
            As = adjacency_matrix.(lattices)
            (Asquare, Acolumnar_dimer, Ahoneycomb) = As

            for (A, lattice) in zip(As, lattices)
                @test tr(A^2)/2 == length(lattice.bonds)
            end


            for n in 2:6
                if any(Ls.<=n)
                    continue
                end
                
                if isodd(n)
                    for A in (Asquare, Acolumnar_dimer, Ahoneycomb)
                        @test tr(A^n) == 0
                    end
                else
                    num_walks = tr(Asquare^n)/prod(Ls)
                    @test num_walks == binomial(n, n÷2)^2
            
                    num_walks = tr(Acolumnar_dimer^n)÷(prod(Ls)*2)
                    @test num_walks == binomial(n, n÷2)^2
                end
            end
        end
    end
end