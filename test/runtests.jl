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
    data = generate_uniform_data(d, n)
end

@testset verbose = true "MondrianForest" begin
    d = 1
    lambda = 4.0
    n_trees = 10
    n = 20
    data = generate_uniform_data(d, n)
    x_eval = ntuple(x -> 0.5, d)
    forest = MondrianForest(lambda, n_trees, x_eval, data["X"], data["Y"])
end
