using MondrianForests
using Test

@testset verbose = true "MondrianCell" begin

    for d in 1:5
        cell = MondrianCell(d)
    end

end

@testset verbose = true "MondrianTree" begin

    cell = MondrianCell(2)
    time = 1.0
    lambda = 5.0
    tree = MondrianTree(cell, time, lambda)

end

