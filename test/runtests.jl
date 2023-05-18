using MondrianForests
using Test
using Distributions

@testset verbose = true "Cells" begin
    for d in 1:5
        cell = MondrianCell(d)
    end
end

@testset verbose = true "Trees" begin
    for d in 1:3
        lambda = 2.0
        tree = MondrianTree(d, lambda)
        show(tree)
    end
end

@testset verbose = true "Data" begin
    d = 2
    n = 10
    data = generate_uniform_data(d, n)
end

@testset verbose = true "Forests" begin
    for d in 1:2
        lambda = 10.0
        n_trees = 1000
        debias_order = 0
        significance_level = 0.05
        n = 2000
        data = generate_uniform_data(d, n)
        x_evals = [ntuple(_ -> x, d) for x in range(0.25, 0.75, length=3)]
        forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                                significance_level, data["X"], data["Y"])
        show(forest)
    end
end

@testset verbose = true "Lifetime selection" begin
    d = 1
    n = 100
    X_dist = X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    eps_dist = Uniform(-sqrt(3), sqrt(3))
    mu = (x -> 3 * x[1]^2)
    sigma2 = (x -> 1/100)
    data = generate_data(n, X_dist, eps_dist, mu, sigma2)
    debias_order = 1
    lambda = select_lifetime_global_polynomial(data["X"], data["Y"], debias_order)
end
