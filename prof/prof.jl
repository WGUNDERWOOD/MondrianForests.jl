using MondrianForests
using Test
using Distributions
using Profile
using ProfileSVG

function prof()
    d = 3
    lambda = 20.0
    n_trees = 2000
    n_data = 10000
    data = generate_uniform_data(d, n_data)
    X_data = data["X"]
    Y_data = data["Y"]
    x_eval = ntuple(i -> 0.5, d)
    forest = MondrianForest(lambda, n_trees, x_eval, X_data, Y_data)
    return nothing
end

prof()
@time prof()
@profile prof()
git_root = strip(read(`git rev-parse --show-toplevel`, String), '\n')
ProfileSVG.save(git_root * "/prof/prof.svg")

