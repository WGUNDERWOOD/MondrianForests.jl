using MondrianForests
using Random
using Distributions

function coverage_experiment(n_reps::Int, lambda::Float64, n_trees::Int, n_data::Int)
    d = 2
    x_eval = ntuple(i -> 0.5, d)
    X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    mu = (x -> x[1]^2 * x[2] + sin(x[1]))
    sigma2 = (x -> x[2]^2 + 1 + x[1])
    eps_dist = Normal(0, 1)
    corrects = fill(false, n_reps)

    Threads.@threads for rep in 1:n_reps
        println(rep)
        data = generate_data(n_data, X_dist, eps_dist, mu, sigma2)
        X_data = data["X"]
        Y_data = data["Y"]
        forest = MondrianForest(lambda, n_trees, x_eval, X_data, Y_data)
        ci_width = 1.96 * sqrt(forest.Sigma_hat) * sqrt(lambda^d / n_data)
        correct = forest.mu_hat - ci_width <= mu(x_eval) <= forest.mu_hat + ci_width
        corrects[rep] = correct
    end
    println("Coverage: ", sum(corrects) / n_reps)
    return nothing
end

n_reps = 100
lambda = 6.0
n_trees = 100
n_data = 1000
coverage_experiment(n_reps, lambda, n_trees, n_data)
