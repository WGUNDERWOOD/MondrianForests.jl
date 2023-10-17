@testset verbose = true "Lifetime polynomial selection" begin
    for debias_order in 0:1
        for d in 1:2
            n = 500
            X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
            eps_dist = Uniform(-sqrt(3), sqrt(3))
            mu = (x -> 3 * x[1]^2)
            sigma2 = (x -> 1 / 100)
            data = MondrianForests.generate_data(n, X_dist, eps_dist, mu, sigma2)
            X_data = data["X"]
            Y_data = data["Y"]
            lambda = select_lifetime_polynomial(X_data, Y_data, debias_order)
        end
    end
end
