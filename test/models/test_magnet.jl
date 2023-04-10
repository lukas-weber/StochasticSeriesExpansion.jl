import StochasticSeriesExpansion.Models as M

@testset "Magnet" begin

    m = M.Magnet(Dict(
        :unitcell => S.UnitCells.square,
        :Lx => 10,
        :Ly => 10,
        :J => 1.2
    ))

end
