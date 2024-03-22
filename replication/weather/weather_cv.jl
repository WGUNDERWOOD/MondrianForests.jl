using CSV
using DataFrames
using PyPlot
using Random
using Colors
using Plots
using Dates
using MondrianForests

# plot setup
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["font.family"] = "serif"
plt.ioff()

function load_data(; limit=nothing)
    # load columns
    file = "replication/weather/weather.csv"
    data = DataFrame(CSV.File(file, limit=limit, missingstring="NA"))
    data = data[:, [:Humidity3pm, :Pressure3pm, :RainTomorrow]]
    dropmissing!(data)
    data.RainTomorrow = replace(data.RainTomorrow, "Yes" => 1, "No" => 0)
    data = Float64.(data)
    # clip outliers
    pres_lims = [985, 1040]
    data = data[pres_lims[1] .<= data.Pressure3pm .<= pres_lims[2], :]
    # jitter
    data.Humidity3pm .+= 1 * (rand(nrow(data)) .- 0.5)
    data.Pressure3pm .+= 0.5 * (rand(nrow(data)) .- 0.5)
    # rescale
    scaling = Dict()
    x_min = minimum(data.Humidity3pm)
    x_max = maximum(data.Humidity3pm)
    y_min = minimum(data.Pressure3pm)
    y_max = maximum(data.Pressure3pm)
    data.Humidity3pm = (data.Humidity3pm .- x_min) ./ (x_max - x_min)
    data.Pressure3pm = (data.Pressure3pm .- y_min) ./ (y_max - y_min)
    Random.seed!(2)
    data = shuffle(data)
    return (data, x_min, x_max, y_min, y_max)
end

function make_evals(n_evals, X, Y)
    eval_ids = sort(shuffle(1:n)[1:n_evals])
    x_evals = X[eval_ids]
    y_evals = Y[eval_ids]
    non_eval_ids = [i for i in 1:n if !(i in eval_ids)]
    X_reduced = X[non_eval_ids]
    Y_reduced = Y[non_eval_ids]
    return (eval_ids, x_evals, y_evals, X_reduced, Y_reduced)
end

# get data and params
(data, x_min, x_max, y_min, y_max) = load_data(limit=nothing)
data = data[1:1000, :]
n = nrow(data)
X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]

n_evals = 200
lambdas = range(0.1, stop=8.0, step=0.1)
n_trees = 400
debias_order = 0

eval_ids = sort(shuffle(1:n)[1:n_evals])
x_evals = X[eval_ids]
y_evals = Y[eval_ids]

# GCV
Random.seed!()
gcvs = zeros(length(lambdas))
mses = zeros(length(lambdas))
d = 2
Threads.@threads for i in 1:length(lambdas)
    lambda = lambdas[i]
    local forest = MondrianForest(lambda, n_trees, x_evals, X, Y)
    mse = sum((y_evals .- forest.mu_hat) .^ 2) / n_evals
    gcv = mse / (1 - lambda^d / n)^2
    gcvs[i] = gcv
    mses[i] = mse
    println("lambda: ", lambda)
    println("gcv: ", gcv)
    println("mse: ", mse)
end

(fig, ax) = plt.subplots(figsize=(4, 4.4))
best_lambda = 5.0
i = [i for i in 1:length(lambdas) if isapprox(lambdas[i], best_lambda, rtol=0.01)][]
plt.plot([best_lambda, best_lambda], [0.0, gcvs[i] - 0.0001], c="#666677",
         linestyle="dashed", lw=1.0)
plt.plot(lambdas, mses, lw=1.0, c="#aa44dd",
         label="Mean squared error")
plt.plot(lambdas, gcvs, lw=1.0, c="#009944",
         label="Generalized cross-validation")
plt.ylim([0.11 - 0.002, 0.19 + 0.002])
plt.yticks(range(0.11, stop=0.19, step=0.01), fontsize=11)
plt.xticks(fontsize=11)
plt.xlabel("Lifetime parameter \$\\lambda\$", fontsize=12)
plt.ylabel("Loss function", fontsize=12)
plt.legend(fontsize=12)
plt.subplots_adjust(left=0.205, right=0.96, top=0.842, bottom=0.140)
plt.savefig("./replication/weather/weather_gcv.png", dpi=350)

sdfsdf

# CIs
limit = nothing
(data, x_min, x_max, y_min, y_max) = load_data(limit=limit)
n = nrow(data)
X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]

n_trees = 400
debias_order = 0
x_evals_original = [(20, 1020), (70, 1000), (80, 990)]
x_evals = [((x[1] - x_min) / (x_max - x_min), (x[2] - y_min) / (y_max - y_min))
           for x in x_evals_original]

forest = MondrianForest(best_lambda, n_trees, x_evals, X, Y, true)

println("Mondrian random forest")
for i in 1:length(x_evals)
    println()
    println("eval: ", x_evals_original[i])
    println("mu_hat: ", forest.mu_hat[i])
    println("CI: ", forest.confidence_band[i])
end

# debiased CIs
n_trees = 200
debias_order = 1
debiased_forest = DebiasedMondrianForest(best_lambda, n_trees, x_evals,
                                         debias_order, X, Y, true)

println("Debiased Mondrian random forest")
for i in 1:length(x_evals)
    println()
    println("eval: ", x_evals_original[i])
    println("mu_hat: ", debiased_forest.mu_hat[i])
    println("CI: ", debiased_forest.confidence_band[i])
end
