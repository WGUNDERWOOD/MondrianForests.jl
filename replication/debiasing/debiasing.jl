using Distributions
using PyPlot
using MondrianForests

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["font.family"] = "serif"
plt.ioff()

d = 1
n = 100000
x_evals = [ntuple(j -> 0.0, d)]
y_evals = [0.0]
n_evals = 1
X_dist = Uniform(-1, 1)
sigma = 0.01
eps_dist = Normal(0, sigma)

rand(X_dist)

X = [ntuple(j -> rand(X_dist), d) for i in 1:n]
Y = [X[i][1]^2 + rand(eps_dist) for i in 1:n]


(fig, ax) = plt.subplots(figsize=(5, 5))
plt.scatter(X, Y)
savefig("replication/debiasing/plot.png", dpi=150)
plt.close()

#lambdas = collect(0.1:0.1:2)
n_trees = 100
#n_subsample = n
#lambda = select_lifetime_gcv(lambdas, n_trees, X, Y, debias_order, n_subsample)
#println(lambda)
for debias_order in [0, 1]
    #lambda = select_lifetime_polynomial(X, Y, debias_order)
    gamma = d + 4*(debias_order + 1)
    println(gamma)
    lambda = n^(-1/gamma)
    println(lambda)
    forest = MondrianForest(lambda, n_trees, x_evals, X, Y)
    mse = forest.mu_hat[]^2
    println(mse)
end

