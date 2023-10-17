using MondrianForests
using PyPlot
using Random
using Distributions

# TODO

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
plt.ioff()
#Random.seed!(314159)

d = 1
X_dist = Uniform(0, 1)
mu = (x -> 10 * (x[1] - 0.5)^2)
sigma2 = (x -> 0.01)
debias_order = 0
eps_dist = Normal(0, 1)
n_data = 1000
n_trees = 1000
significance_level = 0.05

data = generate_data(n_data, X_dist, eps_dist, mu, sigma2)
X_data = data["X"]
Y_data = data["Y"]
lambda_hat = select_lifetime_global_polynomial(X_data, Y_data, debias_order)
#lambda_hat = 10.0
println(lambda_hat)

n_evals = 100
x_evals = [ntuple(_ -> x, 1) for x in range(0, 1, length=n_evals)]

forest = MondrianForest(lambda_hat, n_trees, x_evals, debias_order,
                        significance_level, data["X"], data["Y"])

(fig, ax) = plt.subplots(figsize=(3, 2.1))
plot(x_evals, mu.(x_evals), lw=1, color="k")
plot(x_evals, forest.mu_hat, lw=1, color="r")
scatter(X_data, Y_data, s=1, color="#aaaaaa")

#plt.xticks([0, 0.5, 1])
#plt.yticks([0, 1, 2, 3, 4])
plt.xlabel("\$x_1\$")
plt.ylabel("\$y\$")

git_root = strip(read(`git rev-parse --show-toplevel`, String), '\n')
savefig(git_root * "/replication/plot_mondrian_forest.pdf", bbox_inches="tight")
plt.close("all")
