using Distributions
using MondrianForests

d = 1
n = 200
B = 200
n_reps = 20
X_dist = Uniform(0, 1)
mu = (x -> sum(sin.(pi .* x)))
sigma = 0.3
eps_dist = Normal(0, sigma)
lambda = 10.0
x_evals = [ntuple(j -> 0.5, d)]
J = 0
Sigma_hat = 0.0

t0 = time()
for rep in 1:n_reps
    println(rep)
    X = [ntuple(j -> rand(X_dist), d) for i in 1:n]
    Y = [mu(X[i]) + rand(eps_dist) for i in 1:n]
    forest = DebiasedMondrianForest(lambda, B, x_evals, J, X, Y, true)
    global Sigma_hat += forest.Sigma_hat[] / n_reps
end
println(time() - t0)
println(Sigma_hat)
