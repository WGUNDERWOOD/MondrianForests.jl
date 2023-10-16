using Random
Random.seed!(0)

@testset verbose = true "Forests" begin

    @testset verbose = true "Constant" begin
        n = 50
        lambda = 5.0
        n_trees = 20
        n_evals = 60
        significance_level = 0.05
        estimate_var = false
        for d in 1:3
            X = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> rand(), d) for _ in 1:n_evals]
            Y = [2.0 for _ in 1:n]
            forest = MondrianForest(lambda, n_trees, x_evals, significance_level,
                                    X, Y, estimate_var)
            @suppress show(forest)
            @test all(forest.mu_hat .== 2)
        end
    end

    @testset verbose = true "Linear" begin
        Random.seed!(0)
        n = 2000
        lambda = 20.0
        n_trees = 500
        significance_level = 0.05
        estimate_var = true
        for d in 1:3
            X = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> 0.5, d)]
            Y = sum.(X)
            forest = MondrianForest(lambda, n_trees, x_evals, significance_level,
                                    X, Y, estimate_var)
            ci_lower = [ci[1] for ci in forest.confidence_band]
            ci_upper = [ci[2] for ci in forest.confidence_band]
            @test all(ci_lower .<= forest.mu_hat .<= ci_upper)
            @test isapprox(forest.mu_hat[], 0.5 * d, rtol=0.01)
        end
    end

    @testset verbose = true "Quadratic" begin
        Random.seed!(0)
        n = 2000
        lambda = 20.0
        n_trees = 500
        significance_level = 0.05
        estimate_var = true
        for d in 1:2
            X = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> 0.5, d)]
            Y = sum.(X) .^ 2
            forest = MondrianForest(lambda, n_trees, x_evals, significance_level,
                                    X, Y, estimate_var)
            ci_lower = [ci[1] for ci in forest.confidence_band]
            ci_upper = [ci[2] for ci in forest.confidence_band]
            @test all(ci_lower .<= forest.mu_hat .<= ci_upper)
            @test isapprox(forest.mu_hat[], 0.25 * d ^ 2, rtol=0.01)
        end
    end

end

