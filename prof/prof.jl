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
end

run_mondrian_forest()
@profile run_mondrian_forest()
ProfileSVG.save("prof/prof.svg")

