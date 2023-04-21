using MondrianForests
using Test

@testset verbose = true "MondrianCell" begin
    for d = 1:5
        cell = MondrianCell(d)
    end
end

@testset verbose = true "MondrianTree" begin
    for d = 1:5
        lambda = 10.0
        tree = MondrianTree(d, lambda)
    end

end
