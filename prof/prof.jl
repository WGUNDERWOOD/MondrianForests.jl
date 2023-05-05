using MondrianForests
using Test
using Distributions
using Profile
using ProfileSVG

function prof()
    d = 3
    n_trees = 1000
    n_data = 1000
    debias_order = 2
    for rep in 1:100
        data = generate_uniform_data(d, n_data)
        X_data = data["X"]
        Y_data = data["Y"]
        lambda_hat = select_lifetime_global_polynomial(X_data, Y_data, debias_order)
        x_eval = ntuple(i -> 0.5, d)
        forest = MondrianForest(lambda_hat, n_trees, x_eval, debias_order, data["X"], data["Y"])
    end
    return nothing
end

prof()
@time prof()
@profile prof()
git_root = strip(read(`git rev-parse --show-toplevel`, String), '\n')
ProfileSVG.save(git_root * "/prof/prof.svg")
