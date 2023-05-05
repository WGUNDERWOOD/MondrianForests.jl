using MondrianForests
using Random
using Distributions

function coverage_experiment(n_reps::Int, n_trees::Int, n_data::Int)
    d = 1
    x_eval = ntuple(i -> 0.5, d)
    X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    mu = (x -> 2 + x[1]^2 - sin(10 * x[1]))
    sigma2 = (x -> (x[1] + 2) / 200)
    debias_order = 0
    significance_level = 0.05
    eps_dist = Normal(0, 1)
    corrects = fill(false, n_reps)

    Threads.@threads for rep in 1:n_reps
        println(rep)
        data = generate_data(n_data, X_dist, eps_dist, mu, sigma2)
        X_data = data["X"]
        Y_data = data["Y"]
        lambda_hat = select_lifetime_global_polynomial(X_data, Y_data, debias_order)
        forest = MondrianForest(lambda_hat, n_trees, x_eval, debias_order,
                                significance_level, data["X"], data["Y"])
        correct = forest.confidence_interval[1] <= mu(x_eval) <= forest.confidence_interval[2]
        println(lambda_hat)
        corrects[rep] = correct
    end
    println("Coverage: ", sum(corrects) / n_reps)
    return nothing
end

n_reps = 1000
n_trees = 1000
n_data = 1000
coverage_experiment(n_reps, n_trees, n_data)
