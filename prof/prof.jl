using MondrianForests
using Test
using Distributions
using Profile
using ProfileSVG

function run_mondrian_forest()
    d = 1
    lambda = 20.0
    n_trees = 3000
    n_data = 5000
    data = generate_data(d, n_data)
    X_data = data["X"]
    Y_data = data["Y"]
    x_eval = ntuple(i -> 0.5, d)
    forest = MondrianForest(lambda, n_trees, x_eval, X_data, Y_data)
    #println(forest.X_data)
    #println(forest.cells)
    #println(forest.membership)
    println(forest.mu_hat)
    println(forest.sigma2_hat)
    println(forest.Sigma_hat)
    println()
    return nothing
end

function check_confidence_intervals()
    n_reps = 2000
    d = 1
    lambda = 10.0
    n_trees = 500
    n_data = 3000
    x_eval = ntuple(i -> 0.5, d)
    corrects = fill(false, n_reps)

    Threads.@threads for rep in 1:n_reps
        println(rep)
        data = generate_data(d, n_data)
        X_data = data["X"]
        Y_data = data["Y"]
        forest = MondrianForest(lambda, n_trees, x_eval, X_data, Y_data)
        ci_width = 1.96 * sqrt(forest.Sigma_hat) * sqrt(lambda^d / n_data)
        correct = forest.mu_hat - ci_width <= 0 <= forest.mu_hat + ci_width
        corrects[rep] = correct
    end
    println(sum(corrects) / n_reps)
    return nothing
end

#run_mondrian_forest()
#@profile run_mondrian_forest()
@profile check_confidence_intervals()
#ProfileSVG.save("prof/prof.svg")

