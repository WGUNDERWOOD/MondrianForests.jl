@testset verbose = true "Lifetime selection" begin
    d = 1
    n = 100
    X_dist = X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    eps_dist = Uniform(-sqrt(3), sqrt(3))
    mu = (x -> 3 * x[1]^2)
    sigma2 = (x -> 1 / 100)
    data = generate_data(n, X_dist, eps_dist, mu, sigma2)
    debias_order = 1
    lambda = select_lifetime_global_polynomial(data["X"], data["Y"], debias_order)
end
