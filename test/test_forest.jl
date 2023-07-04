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

