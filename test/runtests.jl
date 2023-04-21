using MondrianForests
using Test
using Distributions

@testset verbose = true "MondrianCell" begin
    for d in 1:5
        cell = MondrianCell(d)
    end
end

@testset verbose = true "MondrianTree" begin
    for d in 1:5
        lambda = 10.0
        tree = MondrianTree(d, lambda)
    end
end

@testset verbose = true "Data" begin
    d = 2
    n = 10
    data = generate_data(d, n)
end

@testset verbose = true "MondrianForest" begin
    d = 2
    lambda = 10.0
    B = 100
    n = 100
    data = generate_data(d, n)
    forest = MondrianForest(d, lambda, B, data["X"], data["Y"])
end
