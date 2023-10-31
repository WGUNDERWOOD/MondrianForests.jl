using Random
Random.seed!(0)

@testset verbose = true "Debiasing" begin
    @testset verbose = true "Constant" begin
        n = 50
        lambda = 5.0
        n_trees = 20
        n_evals = 60
        debias_order = 1
        for d in 1:3
            X_data = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> rand(), d) for _ in 1:n_evals]
            Y_data = [2.0 for _ in 1:n]
            forest = DebiasedMondrianForest(lambda, n_trees, x_evals, debias_order, X_data, Y_data)
            @suppress show(forest)
            @test all(isapprox.(forest.mu_hat, 2, rtol=1e-10))
        end
    end

    @testset verbose = true "Linear" begin
        Random.seed!(0)
        n = 2000
        lambda = 20.0
        n_trees = 500
        debias_order = 1
        estimate_var = true
        for d in 1:2
            X_data = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> 0.5, d)]
            Y_data = sum.(X_data)
            forest = DebiasedMondrianForest(lambda, n_trees, x_evals, debias_order, X_data, Y_data,
                                            estimate_var)
            ci_lower = [ci[1] for ci in forest.confidence_band]
            ci_upper = [ci[2] for ci in forest.confidence_band]
            @test all(ci_lower .<= forest.mu_hat .<= ci_upper)
            @test isapprox(forest.mu_hat[], 0.5 * d, rtol=0.01)
        end
    end

    @testset verbose = true "Quadratic" begin
        Random.seed!(0)
        n = 3000
        lambda = 20.0
        n_trees = 500
        debias_order = 1
        estimate_var = true
        for d in 1:2
            X_data = [ntuple(i -> rand(), d) for _ in 1:n]
            x_evals = [ntuple(i -> 0.5, d)]
            Y_data = sum.(X_data) .^ 2
            forest = DebiasedMondrianForest(lambda, n_trees, x_evals, debias_order, X_data, Y_data,
                                            estimate_var)
            ci_lower = [ci[1] for ci in forest.confidence_band]
            ci_upper = [ci[2] for ci in forest.confidence_band]
            @test all(ci_lower .<= forest.mu_hat .<= ci_upper)
            @test isapprox(forest.mu_hat[], 0.25 * d^2, rtol=0.01)
        end
    end
end
