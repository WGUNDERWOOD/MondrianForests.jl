using Distributions

"""
    generate_data(n::Int, X_dist::Distribution, eps_dist::Distribution,
                  mu::Function, sigma2::Function)

Generate sample data for Mondrian forest estimation.

Draws `n` independent samples from \$Y = \\mu(X) + \\sigma(X) \\varepsilon\$,
with \$X \\sim\$ `X_dist` and \$\\varepsilon \\sim\$ `eps_dist`.
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
    generate_uniform_data_uniform_errors(d::Int, n::Int)

Generate uniform sample data with uniform errors for Mondrian forest estimation.

Draws `n` independent samples from \$Y = \\varepsilon\$,
with \$X \\sim \\mathcal{U}[0, 1]\$ and
\$\\varepsilon \\sim \\mathcal{U}\\big[-\\sqrt 3, \\sqrt 3\\big]\$.
"""
function generate_uniform_data_uniform_errors(d::Int, n::Int)
    X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    eps_dist = Uniform(-sqrt(3), sqrt(3))
    mu = (x -> 0.0)
    sigma2 = (x -> 1.0)
    return generate_data(n, X_dist, eps_dist, mu, sigma2)
end

"""
    generate_uniform_data_normal_errors(d::Int, n::Int)

Generate uniform sample data with normal errors for Mondrian forest estimation.

Draws `n` independent samples from \$Y = \\varepsilon\$,
with \$X \\sim \\mathcal{U}[0, 1]\$ and \$\\varepsilon \\sim \\mathcal{N}(0, 1)\$.
"""
function generate_uniform_data_normal_errors(d::Int, n::Int)
    X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    eps_dist = Normal(0, 1)
    mu = (x -> 0.0)
    sigma2 = (x -> 1.0)
    return generate_data(n, X_dist, eps_dist, mu, sigma2)
end
