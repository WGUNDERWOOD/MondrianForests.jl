using MondrianForests
using PyPlot
using Random

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
plt.ioff()

d = 1
x_eval = ntuple(i -> 0.5, d)
X_dist = Beta(2, 3)
mu = (x -> x[1]^2 + sin(x[1]))
sigma2 = (x -> (x[1] + 1) / 10)
debias_order = 1
eps_dist = Normal(0, 1)

data = generate_data(n_data, X_dist, eps_dist, mu, sigma2)
X_data = data["X"]
Y_data = data["Y"]
lambda_hat = select_lifetime_global_polynomial(X_data, Y_data, debias_order)
forest = MondrianForest(lambda_hat, n_trees, x_eval, debias_order, data["X"], data["Y"])

lambda = lambdas[i]
tree = MondrianTree(d, lambda)

#(fig, ax) = plot_mondrian_tree(tree)
#git_root = strip(read(`git rev-parse --show-toplevel`, String), '\n')
#savefig(git_root * "/replication/plot_mondrian_process_$i.pgf", bbox_inches="tight")
#plt.close("all")
