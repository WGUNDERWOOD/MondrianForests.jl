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

@testset verbose = true "Estimation" begin
    d = 1
    lambda = 10.0
    n_trees = 500
    n = 2000
    data = generate_data(d, n)
    forest = MondrianForest(d, lambda, n_trees)
    x = [0.5 for _ in 1:d]
    MondrianForests.fit(forest, x, data["X"], data["Y"])
    println(forest)
    println(0.409^d)
end
