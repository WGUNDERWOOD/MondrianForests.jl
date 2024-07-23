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
    gcv
end

mutable struct Experiment
    # parameters
    J_estimator::Int
    J_lifetime::Int
    lambda_target::LambdaTarget
    lambda_method::LambdaMethod
    lambda_multipliers::Vector{Float64}
    lambda_candidates::Vector{Float64}
    n_subsample::Int

    # distribution
    d::Int
    n::Int
    B_estimator::Int
    B_lifetime::Int
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
    sd_theory::Float64
    bias_theory::Float64
end

function Experiment(
        J_estimator::Int,
        lambda_target::LambdaTarget,
        lambda_method::LambdaMethod,
        lambda_multipliers::Vector{Float64},
        lambda_candidates::Vector{Float64},
        n_subsample::Int,
        d::Int,
        n::Int,
        B_estimator::Int,
        B_lifetime::Int,
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
               lambda_candidates,
               n_subsample,
               d,
               n,
               B_estimator,
               B_lifetime,
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
               NaN,
               NaN
              )
end

function run_first_block()
    J_estimator = 0
    lambda_target = rmse::LambdaTarget
    lambda_methods = instances(LambdaMethod)
    lambda_multipliers = [0.9, 1.0]
    lambda_candidates = [4.0, 5.0]
    n_subsample = 20
    d = 1
    n = 1000
    B_estimator = 200
    B_lifetime = 200
    x_evals = [ntuple(j -> 0.5, d)]
    X_dist = Uniform(0, 1)
    mu = (x -> sum(sin.(x)))
    sigma = 0.1
    eps_dist = Normal(0, sigma)
    experiment = Experiment(J_estimator, lambda_target,
                            lambda_methods[1], lambda_multipliers,
                            lambda_candidates, n_subsample,
                            d, n, B_estimator,
                            B_lifetime, x_evals, X_dist, mu, eps_dist)
    run(experiment)
end

function get_J_lifetime(experiment)
    if experiment.lambda_target == rmse::LambdaTarget
        return experiment.J_estimator
    elseif experiment.lambda_target == undersmooth::LambdaTarget
        return experiment.J_estimator - 1
    end
end

function get_theory(experiment::Experiment)
    n = experiment.n
    d = experiment.d
    lambda = experiment.lambda
    x_evals = experiment.x_evals
    sigma2 = var(experiment.eps_dist)
    if experiment.J_estimator == 0
        C0 = (4 - 4 * log(2)) / 3
        experiment.sd_theory = sqrt(lambda^d * sigma2 * C0^d / n)
        # TODO calculate sin bias
        experiment.bias_theory = - sum(sin.(x_evals[])) / (2 * lambda^2)
    elseif experiment.J_estimator == 1
        C1 = (4/3 - 4*log(2)/3)
        C2 = (2 - 2*log(2))
        C3 = (5/3 - log(5/2) - 3*log(5/3)/2)
        C_all = 16/5 * C1^d + 81/25 * C2^d - 72/5 * C3^d
        experiment.sd_theory = sqrt(lambda^d * sigma2 * C_all / n)
        # TODO calculate sin bias
        experiment.bias_theory = d / lambda^2
    end
end

function select_lifetime(X, Y, experiment)
    J = experiment.J_lifetime
    lambda_candidates = experiment.lambda_candidates
    B = experiment.B_lifetime
    n_subsample = experiment.n_subsample
    if experiment.lambda_method == optimal::LambdaMethod
        # TODO
        return 10.0
    elseif experiment.lambda_method == polynomial::LambdaMethod
        return select_lifetime_polynomial(X, Y, J)
    elseif experiment.lambda_method == gcv::LambdaMethod
        return select_lifetime_gcv(lambdas, n_trees, X, Y, J)
    end
end

function run(experiment::Experiment)
    n_rep = 200
    n = experiment.n
    d = experiment.d
    x_evals = experiment.x_evals
    mu = experiment.mu
    experiment.J_lifetime = get_J_lifetime(experiment)
    mse = 0.0
    bias = 0.0
    coverage = 0.0
    average_width = 0.0
    average_lambda = 0.0
    for rep in 1:n_rep 
        println(rep)
        X = [ntuple(j -> rand(experiment.X_dist), d) for i in 1:n]
        Y = [mu(X[i]) + rand(experiment.eps_dist) for i in 1:n]
        lambda = select_lifetime(X, Y, experiment)
        forest = DebiasedMondrianForest(lambda, experiment.B_estimator,
                                        x_evals,
                                        experiment.J_estimator, X, Y, true)
        #println(forest.debias_scaling)
        #println(forest.debias_coeffs)
        ci = forest.confidence_band
        mse += (forest.mu_hat[] - mu(x_evals[]))^2 / n_rep
        bias += (forest.mu_hat[] - mu(x_evals[])) / n_rep
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
    get_theory(experiment)
    #show(experiment)

    for f in fieldnames(Experiment)
        v = getfield(experiment, f)
        println("$f: $v")
    end

end

run_first_block()

# params
#d = 1
#n = 50
#x_evals = [ntuple(j -> 0.0, d)]
#y_evals = [0.0]
#n_evals = 1
#X_dist = Uniform(-1, 1)
#sigma = 0.001
#eps_dist = Normal(0, sigma)
#X = [ntuple(j -> rand(X_dist), d) for i in 1:n]
#Y = [X[i][1]^2 + rand(eps_dist) for i in 1:n]

# plot data
#(fig, ax) = plt.subplots(figsize=(5, 5))
#plt.scatter(X, Y)
#savefig("replication/debiasing/plot.png", dpi=150)
#plt.close()


# run experiment
#lambdas = collect(1:0.1:4)
#n_trees = 100
#n_reps = 50
#n_subsample = n
#lambda = select_lifetime_gcv(lambdas, n_trees, X, Y, debias_order, n_subsample)
#println(lambda)
