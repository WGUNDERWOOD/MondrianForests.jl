using MondrianForests
using Random
using Distributions

function coverage_experiment(n_reps::Int, n_trees::Int, n_data::Int)
    d = 1
    x_evals = [ntuple(i -> 0.5, d)]
    X_dist = product_distribution([Uniform(0, 1) for _ in 1:d])
    #mu = (x -> 2 + x[1]^2 - sin(10 * x[1]))
    mu = (x -> 10 * (x[1] - 0.5)^2)
    #mu = (x -> 1)
    #sigma2 = (x -> (x[1] + 2) / 200)
    sigma2 = (x -> 0.01)
    debias_order = 0
    significance_level = 0.05
    eps_dist = Normal(0, 1)
    corrects0 = fill(false, n_reps)
    corrects1 = fill(false, n_reps)

    Threads.@threads for rep in 1:n_reps
        println(rep)
        data = generate_data(n_data, X_dist, eps_dist, mu, sigma2)
        X_data = data["X"]
        Y_data = data["Y"]
        lambda_hat = select_lifetime_global_polynomial(X_data, Y_data, debias_order)
        #lambda_hat = 20.0
        forest0 = MondrianForest(lambda_hat, n_trees, x_evals, debias_order,
                                 significance_level, data["X"], data["Y"])
        forest1 = MondrianForest(lambda_hat, n_trees, x_evals, debias_order + 1,
                                 significance_level, data["X"], data["Y"])

        correct0 = all(forest0.confidence_band[s][1] <= mu(x_evals[s]) <=
                       forest0.confidence_band[s][2] for s in 1:(forest0.n_evals))
        correct1 = all(forest1.confidence_band[s][1] <= mu(x_evals[s]) <=
                       forest1.confidence_band[s][2] for s in 1:(forest1.n_evals))
        #
        #println(forest0.sigma2_hat)
        #println(forest1.sigma2_hat)
        corrects0[rep] = correct0
        corrects1[rep] = correct1
    end
    println("Coverage 0: ", sum(corrects0) / n_reps)
    println("Coverage 1: ", sum(corrects1) / n_reps)
    return nothing
end

n_reps = 1000
n_trees = 500
n_data = 500
coverage_experiment(n_reps, n_trees, n_data)
