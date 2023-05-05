using Distributions

"""
Generate sample data for Mondrian forest estimation.

Draws `n` independent samples from `Y = mu(X) + sigma(X) eps`,
with `X ~ X_dist` and `eps ~ eps_dist`
"""
function generate_data(n::Int, X_dist::Distribution, eps_dist::Distribution,
                       mu::Function, sigma2::Function)
    d = length(X_dist)
    X = [rand(X_dist) for _ in 1:n]
    X = [ntuple(j -> X[i][j], Val(d)) for i in 1:n]
    eps = [rand(eps_dist) for _ in 1:n]
    Y = mu.(X) + sqrt.(sigma2.(X)) .* eps
    return Dict("X" => X, "Y" => Y)
end

"""
Generate uniform sample data for Mondrian forest estimation.

Draws `n` independent samples from `Y = eps`,
with `X ~ U[0, 1]` and `eps ~ U[-sqrt(3), sqrt(3)]`
"""
function generate_uniform_data(d::Int, n::Int)
    X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    eps_dist = Uniform(-sqrt(3), sqrt(3))
    mu = (x -> 0.0)
    sigma2 = (x -> 1.0)
    return generate_data(n, X_dist, eps_dist, mu, sigma2)
end
