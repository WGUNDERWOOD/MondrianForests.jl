@testset verbose = true "Forests" begin

    @testset verbose = true "Constant" begin
        n = 50
        lambda = 5.0
        n_trees = 20
        n_evals = 60
        debias_order = 0
        significance_level = 0.05
        for d in 1:3
            X = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> rand(), d) for _ in 1:n_evals]
            Y = [2.0 for _ in 1:n]
            forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                                    significance_level, X, Y)
            @suppress show(forest)
            @test all(forest.mu_hat .== 2)
        end
    end

    @testset verbose = true "Linear" begin
        Random.seed!(0)
        n = 5000
        lambda = 20.0
        n_trees = 100
        debias_order = 0
        significance_level = 0.05
        for d in 1:3
            X = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> 0.5, d)]
            Y = sum.(X)
            forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                                    significance_level, X, Y)
            @test isapprox(forest.mu_hat[], 0.5 * d, rtol=0.01)
        end
    end

    @testset verbose = true "Quadratic" begin
        Random.seed!(0)
        n = 5000
        lambda = 20.0
        n_trees = 1000
        debias_order = 0
        significance_level = 0.05
        for d in 1:2
            X = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> 0.5, d)]
            Y = sum.(X) .^ 2
            forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                                    significance_level, X, Y)
            @test isapprox(forest.mu_hat[], 0.25 * d ^ 2, rtol=0.01)
        end
    end

    @testset verbose = true "Quadratic debiased" begin
        # TODO this is failing with the debiasing
        Random.seed!(0)
        n = 5000
        lambda = 20.0
        n_trees = 1000
        debias_order = 1
        significance_level = 0.05
        for d in 1:2
            X = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> 0.5, d)]
            Y = sum.(X) .^ 2
            forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                                    significance_level, X, Y)
            @test isapprox(forest.mu_hat[], 0.25 * d ^ 2, rtol=0.01)
        end
    end

end

