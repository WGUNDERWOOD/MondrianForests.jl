using Distributions
using MondrianForests
using DataFrames
using CSV
using Random

@enum LifetimeMethod begin
    opt
    pol
end

mutable struct Experiment
    # parameters
    J_estimator::Int
    J_lifetime::Int
    lifetime_method::LifetimeMethod
    lifetime_multiplier::Float64

    # distribution
    d::Int
    n::Int
    B::Int
    x_evals
    X_dist::Distribution
    mu::Function
    eps_dist::Distribution

    # data
    # TODO remove these to reduce memory use
    X
    Y
    rep::Int

    # outputs
    mu_hat::Float64
    sd_hat::Float64
    sigma2_hat::Float64
    coverage::Bool
    width::Float64
    lambda::Float64
    rmse_theory::Float64
    sd_theory::Float64
    bias_theory::Float64
end

function Experiment(J_estimator::Int, J_lifetime::Int, lifetime_method::LifetimeMethod,
        lifetime_multiplier::Float64, d::Int, n::Int, B::Int, x_evals,
        X_dist::Distribution, mu::Function, eps_dist::Distribution, X, Y, rep)
    Experiment(J_estimator, J_lifetime, lifetime_method, lifetime_multiplier, d, n, B,
               x_evals, X_dist, mu, eps_dist, X, Y, rep, NaN, NaN, NaN, false, NaN,
               NaN, NaN, NaN, NaN)
end

function run_all()
    # tables format is (d, n, B)
    tables = [
              (1, 1000, 800), # good
              (2, 1000, 800), # good
              #(1, 1000, 500), # medium test
              #(2, 1000, 200), # medium test
              #(1, 100, 100), # small test
              #(2, 10, 10), # small test
             ]
    n_reps = 3000
    lifetime_methods = [opt::LifetimeMethod, pol::LifetimeMethod]
    #lifetime_methods = [opt::LifetimeMethod]
    lifetime_multipliers = [0.8, 0.9, 1.0, 1.1, 1.2]
    #lifetime_multipliers = [1.0]
    X_dist = Uniform(0, 1)
    mu = (x -> sum(sin.(pi .* x)))
    sigma = 0.3
    eps_dist = Normal(0, sigma)
    J_blocks = [(0, 0), (1, 1), (1, 0)]
    #J_blocks = [(1, 0)]
    experiments = []

    for (d, n, B) in tables
        x_evals = [ntuple(j -> 0.5, d)]
        for rep in 1:n_reps
            X = [ntuple(j -> rand(X_dist), d) for i in 1:n]
            Y = [mu(X[i]) + rand(eps_dist) for i in 1:n]
            for (J_estimator, J_lifetime) in J_blocks
                for lifetime_method in lifetime_methods
                    if lifetime_method == opt::LifetimeMethod
                        lifetime_mults = lifetime_multipliers
                    else
                        lifetime_mults = [1.0]
                    end
                    for lifetime_multiplier in lifetime_mults
                        experiment = Experiment(J_estimator, J_lifetime, lifetime_method,
                                                lifetime_multiplier, d, n, B,
                                                x_evals, X_dist, mu, eps_dist, X, Y, rep)
                        push!(experiments, experiment)
                    end
                end
            end
        end
    end

    shuffle!(experiments)
    count = 1
    t0 = time()
    n_exp = length(experiments)
    Threads.@threads for experiment in experiments
        f = "d = $(experiment.d), n = $(experiment.n), B = $(experiment.B), "
        f *= "Je = $(experiment.J_estimator), Jl = $(experiment.J_lifetime), "
        f *= "rep = $(experiment.rep)"
        t1 = time() - t0
        rate = count / t1
        t_left = (n_exp - count) / rate
        println(round(t_left, digits=0), "s left, ",
                round(t_left / 60, digits=2), "min left")
        println(f)
        println("$count / $n_exp")
        count += 1
        run(experiment)
    end

    # summarize the results of each experiment
    results = []
    for (d, n, B) in tables
        for (J_estimator, J_lifetime) in J_blocks
            for lifetime_method in instances(LifetimeMethod)
                if lifetime_method == opt::LifetimeMethod
                    lifetime_multipliers = [0.8, 0.9, 1.0, 1.1, 1.2]
                else
                    lifetime_multipliers = [1.0]
                end
                for lifetime_multiplier in lifetime_multipliers
                    experiments_small = [e for e in experiments if
                                         (e.d, e.n, e.B, e.J_estimator, e.J_lifetime,
                                          e.lifetime_method, e.lifetime_multiplier)
                                         == (d, n, B, J_estimator, J_lifetime, lifetime_method,
                                             lifetime_multiplier)]
                    n_small = length(experiments_small)
                    if n_small > 0
                        result = Dict(
                                      "d" => d,
                                      "n" => n,
                                      "B" => B,
                                      "J_estimator" => J_estimator,
                                      "J_lifetime" => J_lifetime,
                                      "lifetime_method" => lifetime_method,
                                      "lifetime_multiplier" => lifetime_multiplier,
                                      "lambda" => sum(e.lambda for e in experiments_small) / n_small,
                                      "rmse" => sqrt(sum((e.mu_hat - e.mu(e.x_evals[]))^2 for e in experiments_small) / n_small),
                                      "bias" => sum(e.mu_hat - e.mu(e.x_evals[]) for e in experiments_small) / n_small,
                                      "sd_hat" => sum(e.sd_hat for e in experiments_small) / n_small,
                                      "sigma2_hat" => sum(e.sigma2_hat for e in experiments_small) / n_small,
                                      "bias_theory" => sum(e.bias_theory for e in experiments_small) / n_small,
                                      "sd_theory" => sum(e.sd_theory for e in experiments_small) / n_small,
                                      "rmse_theory" => sum(e.rmse_theory for e in experiments_small) / n_small,
                                      "coverage" => sum(e.coverage for e in experiments_small) / n_small,
                                      "average_width" => sum(e.width for e in experiments_small) / n_small,
                                     )
                        result["sd"] = sqrt(result["rmse"]^2 - result["bias"]^2)
                        result["bias_over_sd"] = abs(result["bias"]) / result["sd"]
                        push!(results, result)
                    end
                end
            end
        end
    end

    df = DataFrame(results)
    #display(df)
    CSV.write("./replication/debiasing/results.csv", df)
end


function get_theory(experiment::Experiment)
    n = experiment.n
    d = experiment.d
    lambda = experiment.lambda
    x_evals = experiment.x_evals
    sigma2 = var(experiment.eps_dist)
    if experiment.J_estimator == 0
        experiment.sd_theory = sqrt(lambda^d * sigma2 * 0.4091^d / n)
        experiment.bias_theory = - pi^2 * d / (2 * lambda^2)
    elseif experiment.J_estimator == 1
        C = 3.2 * 0.4091^d - 2.88 * 0.4932^d + 3.24 * 0.6137^d
        experiment.sd_theory = sqrt(lambda^d * sigma2 * C / n)
        experiment.bias_theory = -4 * pi^4 * d / (27 * lambda^4)
    end
    experiment.rmse_theory = sqrt(experiment.bias_theory^2 + experiment.sd_theory^2)
end

function select_lifetime(X, Y, x_eval, experiment)
    d = experiment.d
    n = experiment.n
    sigma2 = var(experiment.eps_dist)
    J_lifetime = experiment.J_lifetime
    if experiment.lifetime_method == opt::LifetimeMethod
        if J_lifetime == 0
            numerator = d * pi^4 * n
            denominator = sigma2 * 0.4091^d
            return (numerator / denominator)^(1 / (4+d))
        elseif J_lifetime == 1
            C = 3.2 * 0.4091^d - 2.88 * 0.4932^d + 3.24 * 0.6137^d
            numerator = 128 * d * pi^8 * n
            denominator = 27^2 * sigma2 * C
            return (numerator / denominator)^(1 / (8+d))
        end
    elseif experiment.lifetime_method == pol::LifetimeMethod
        return select_lifetime_polynomial_amse(X, Y, x_eval, J_lifetime)
    end
end

function run(experiment::Experiment)
    n = experiment.n
    d = experiment.d
    B = experiment.B
    J_estimator = experiment.J_estimator
    x_evals = experiment.x_evals
    X = experiment.X
    Y = experiment.Y
    mu = experiment.mu
    lifetime_multiplier = experiment.lifetime_multiplier
    lambda = select_lifetime(X, Y, x_evals[], experiment) * lifetime_multiplier
    forest = DebiasedMondrianForest(lambda, B, x_evals, J_estimator, X, Y, true)
    experiment.mu_hat = forest.mu_hat[]
    ci = forest.confidence_band
    experiment.sd_hat = sqrt(forest.Sigma_hat[] * lambda^d / n)
    experiment.sigma2_hat = forest.sigma2_hat[]
    experiment.coverage = (ci[][1] <= mu(x_evals[]) <= ci[][2])
    experiment.width = ci[][2] - ci[][1]
    experiment.lambda = lambda
    get_theory(experiment)
end

run_all()

# TODO use medians?
