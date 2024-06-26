using MondrianForests

# sample a two-dimensional Mondrian tree
d = 2
lambda = 2.0
tree = MondrianTree(d, lambda)
println()
show(tree)
println()

# generate some data
# covariates X_data are two-dimensional
# response Y_data is one-dimensional
# true regression function is zero
n_data = 100
data = MondrianForests.generate_uniform_data_uniform_errors(d, n_data)
X_data = data["X"]
Y_data = data["Y"]
println("covariates: ")
display(X_data[1:5])
println("\nresponses: ")
display(Y_data[1:5])

# select a lifetime parameter
# with generalized cross-validation
n_trees = 50
n_subsample = 30
debias_order = 0
lambdas = collect(range(0.5, 10.0, step=0.5))
lambda = select_lifetime_gcv(lambdas, n_trees, X_data, Y_data, debias_order, n_subsample)
println("\nlambda chosen by GCV: ", lambda)

# fit and evaluate a Mondrian random forest
x_evals = [(0.5, 0.5), (0.2, 0.8)]
estimate_var = true
forest = MondrianForest(lambda, n_trees, x_evals, X_data, Y_data, estimate_var)
println("\nestimated regression function:")
display(forest.mu_hat)
println("\nestimated estimator variance:")
display(forest.Sigma_hat)
println("\nestimated confidence band:")
display(forest.confidence_band)

# fit and evaluate a debiased Mondrian random forest
debiased_forest = DebiasedMondrianForest(lambda, n_trees, x_evals, debias_order,
                                         X_data, Y_data, estimate_var)
println("\ndebiased estimated regression function:")
display(debiased_forest.mu_hat)
println("\ndebiased estimated estimator variance:")
display(debiased_forest.Sigma_hat)
println("\ndebiased estimated confidence band:")
display(debiased_forest.confidence_band)
