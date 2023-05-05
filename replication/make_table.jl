using MondrianForests
using Random
using Distributions

function coverage_experiment(n_reps::Int, n_trees::Int, n_data::Int)
    d = 2
    x_eval = ntuple(i -> 0.5, d)
    X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    mu = (x -> x[1]^2 * x[2] + sin(x[1]))
    sigma2 = (x -> x[2]^2 + 1 + x[1])
    debias_order = 1
    eps_dist = Normal(0, 1)
    corrects = fill(false, n_reps)

    for rep in 1:n_reps
        println(rep)
        data = generate_data(n_data, X_dist, eps_dist, mu, sigma2)
        X_data = data["X"]
        Y_data = data["Y"]
        lambda_hat = select_lifetime_global_polynomial(X_data, Y_data, debias_order)
        forest = MondrianForest(lambda_hat, n_trees, x_eval, debias_order, data["X"], data["Y"])
        ci_width = 1.96 * sqrt(forest.Sigma_hat) * sqrt(lambda_hat^d / n_data)
        correct = forest.mu_hat - ci_width <= mu(x_eval) <= forest.mu_hat + ci_width
        corrects[rep] = correct
    end
    println("Coverage: ", sum(corrects) / n_reps)
    return nothing
end

n_reps = 1000
n_trees = 1000
n_data = 1000
coverage_experiment(n_reps, n_trees, n_data)
