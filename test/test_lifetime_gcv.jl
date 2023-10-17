@testset verbose = true "Lifetime GCV selection" begin
    for debias_order in 0:1
        for d in 1:2
            n = 50
            X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
            eps_dist = Uniform(-sqrt(3), sqrt(3))
            mu = (x -> 3 * x[1]^2)
            sigma2 = (x -> 1.0)
            data = MondrianForests.generate_data(n, X_dist, eps_dist, mu, sigma2)
            X_data = data["X"]
            Y_data = data["Y"]


            n_trees = 100
            n_subsample = 20
            lambdas = collect(range(1.0, 20.1, step=1))
            lambda = select_lifetime_gcv(lambdas, n_trees, n_subsample, X_data, Y_data, debias_order)
            println(lambdas)
            println(lambda)
        end
    end
end
