using Distributions
using PyPlot
using MondrianForests

# plot setup
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["font.family"] = "serif"
plt.ioff()

# params
d = 1
n = 50
x_evals = [ntuple(j -> 0.0, d)]
y_evals = [0.0]
n_evals = 1
X_dist = Uniform(-1, 1)
sigma = 0.001
eps_dist = Normal(0, sigma)


# plot data
#(fig, ax) = plt.subplots(figsize=(5, 5))
#plt.scatter(X, Y)
#savefig("replication/debiasing/plot.png", dpi=150)
#plt.close()

# run experiment
lambdas = collect(1:0.1:4)
n_trees = 100
n_reps = 50
#n_subsample = n
#lambda = select_lifetime_gcv(lambdas, n_trees, X, Y, debias_order, n_subsample)
#println(lambda)
for debias_order in [0, 1]
    mses = Float64[]
    for rep in 1:n_reps
        X = [ntuple(j -> rand(X_dist), d) for i in 1:n]
        Y = [X[i][1]^2 + rand(eps_dist) for i in 1:n]
        lambda = select_lifetime_polynomial(X, Y, debias_order)
        #lambda = select_lifetime_gcv(lambdas, n_trees, X, Y, debias_order, n)
        println(lambda)
        forest = MondrianForest(lambda, n_trees, x_evals, X, Y)
        mse = forest.mu_hat[]^2
        push!(mses, mse)
        #println(lambda)
    end
    mse_mean = sum(mses) / length(mses)
    #println(mses)
    println(mse_mean)
end
