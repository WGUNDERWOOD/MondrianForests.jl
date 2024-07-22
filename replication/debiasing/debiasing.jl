using Distributions
using PyPlot
using MondrianForests

# plot setup
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["font.family"] = "serif"
plt.ioff()

@enum LambdaTarget begin
    rmse
    undersmooth
end

@enum LambdaMethod begin
    optimal
    polynomial
    loocv
    gcv
end

mutable struct Experiment
    # parameters
    J_estimator::Int
    J_lifetime::Int
    lambda_target::LambdaTarget
    lambda_method::LambdaMethod
    lambda_multipliers::Vector{Float64}

    # distribution
    d::Int
    n::Int
    B::Int
    x_evals
    X_dist::Distribution
    mu::Function
    eps_dist::Distribution

    # outputs
    rmse::Float64
    bias::Float64
    sd::Float64
    bias_over_sd::Float64
    coverage::Float64
    average_width::Float64
    lambda::Float64
end

function Experiment(
        J_estimator::Int,
        lambda_target::LambdaTarget,
        lambda_method::LambdaMethod,
        lambda_multipliers::Vector{Float64},
        d::Int,
        n::Int,
        B::Int,
        x_evals,
        X_dist::Distribution,
        mu::Function,
        eps_dist::Distribution,
    )
    Experiment(
               J_estimator,
               0,
               lambda_target,
               lambda_method,
               lambda_multipliers,
               d,
               n,
               B,
               x_evals,
               X_dist,
               mu,
               eps_dist,
               NaN,
               NaN,
               NaN,
               NaN,
               NaN,
               NaN,
               NaN,
              )
end

function run_first_block()
    J_estimator = 0
    lambda_target = rmse::LambdaTarget
    lambda_methods = instances(LambdaMethod)
    lambda_multipliers = [0.9, 1.0]
    d = 1
    n = 20
    B = 20
    x_evals = [ntuple(j -> 0.0, d)]
    X_dist = Uniform(-1, 1)
    mu = (x -> x[1]^2)
    sigma = 0.01
    eps_dist = Normal(0, sigma)
    experiment = Experiment(J_estimator, lambda_target,
                            lambda_methods[1], lambda_multipliers,
                            d, n, B, x_evals, X_dist, mu, eps_dist)
    run(experiment)
end

function get_J_lifetime(experiment)
    if experiment.lambda_target == rmse::LambdaTarget
        return experiment.J_estimator
    elseif experiment.lambda_target == undersmooth::LambdaTarget
        return experiment.J_estimator - 1
    end
end

function run(experiment::Experiment)
    n_rep = 10
    n = experiment.n
    B = experiment.B
    d = experiment.d
    x_evals = experiment.x_evals
    experiment.J_lifetime = get_J_lifetime(experiment)
    mse = 0.0
    bias = 0.0
    coverage = 0.0
    average_width = 0.0
    average_lambda = 0.0
    for rep in 1:n_rep 
        println(rep)
        X = [ntuple(j -> rand(experiment.X_dist), d) for i in 1:n]
        Y = [X[i][1]^2 + rand(experiment.eps_dist) for i in 1:n]
        lambda = select_lifetime_polynomial(X, Y, experiment.J_lifetime)
        forest = DebiasedMondrianForest(lambda, experiment.B, x_evals,
                                        experiment.J_estimator, X, Y, true)
        ci = forest.confidence_band
        mse += forest.mu_hat[]^2 / n_rep
        bias += forest.mu_hat[] / n_rep
        coverage += (ci[][1] <= 0 <= ci[][2]) / n_rep
        average_width += (ci[][2] - ci[][1]) / n_rep
        average_lambda += lambda / n_rep
    end
    experiment.rmse = sqrt(mse)
    experiment.bias = bias
    experiment.sd = sqrt(mse - bias^2)
    experiment.bias_over_sd = abs(bias) / experiment.sd
    experiment.coverage = coverage
    experiment.average_width = average_width
    experiment.lambda = average_lambda
    show(experiment)
end

run_first_block()

# params
d = 1
n = 50
x_evals = [ntuple(j -> 0.0, d)]
y_evals = [0.0]
n_evals = 1
X_dist = Uniform(-1, 1)
sigma = 0.001
eps_dist = Normal(0, sigma)
X = [ntuple(j -> rand(X_dist), d) for i in 1:n]
Y = [X[i][1]^2 + rand(eps_dist) for i in 1:n]

# plot data
(fig, ax) = plt.subplots(figsize=(5, 5))
plt.scatter(X, Y)
savefig("replication/debiasing/plot.png", dpi=150)
plt.close()

asdfsdf

# run experiment
lambdas = collect(1:0.1:4)
n_trees = 100
n_reps = 50
#n_subsample = n
#lambda = select_lifetime_gcv(lambdas, n_trees, X, Y, debias_order, n_subsample)
#println(lambda)
