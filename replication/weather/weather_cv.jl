using Revise
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
rcParams["text.latex.preamble"]="\\usepackage[sfdefault,light]{FiraSans}"
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
    Random.seed!(1)
    data = shuffle(data)
    #println("number of samples ", nrow(data))
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

function get_mse_fold(n_evals, X, Y, lambda, n_trees, debias_order)
    (eval_ids, x_evals, y_evals, X_reduced, Y_reduced) = make_evals(n_evals, X, Y)
    forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                            0.05, X_reduced, Y_reduced, false)
    mse_fold = sum((y_evals .- forest.mu_hat).^2) / n_evals
    return mse_fold
end

function get_mse(n_evals, X, Y, lambda, n_trees, debias_order, n_folds)
    mse = 0.0
    for fold in 1:n_folds
        println("fold: ", fold)
        mse_fold = get_mse_fold(n_evals, X, Y, lambda, n_trees, debias_order)
        mse += mse_fold / n_folds
    end
    return mse
end

# get data and params
(data, x_min, x_max, y_min, y_max) = load_data(limit=nothing)
data = data[1:1000, :]
n = nrow(data)
X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]

n_evals = 200
lambdas = range(0.1, stop=8.0, step=0.1)
n_trees = 200
debias_order = 0

eval_ids = sort(shuffle(1:n)[1:n_evals])
x_evals = X[eval_ids]
y_evals = Y[eval_ids]

# GCV
Random.seed!()
gcvs = Float64[]
mses = Float64[]
d = 2
for lambda in lambdas
    println("lambda: ", lambda)
    forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                            0.05, X, Y, false, false)
    #gcv_dof = forest.gcv_dof
    mse = sum((y_evals .- forest.mu_hat).^2) / n_evals
    #gcv = mse / ((1 - gcv_dof / n)^2)
    gcv = mse / (1 - lambda^d / n)^2
    push!(gcvs, gcv)
    push!(mses, mse)
    #println("gcv dof: ", gcv_dof)
    println("gcv: ", gcv)
    println("mse: ", mse)
end

(fig, ax) = plt.subplots(figsize=(3.5, 3.5))
best_lambda = 5.0
i = [i for i in 1:length(lambdas) if lambdas[i] == best_lambda][]
plt.plot([best_lambda, best_lambda], [0.0, gcvs[i] - 0.0001], c="#666677",
         linestyle="dashed", lw=1.0)
plt.plot(lambdas, mses, lw=1.0, c="#aa44dd",
         label="Mean squared error")
plt.plot(lambdas, gcvs, lw=1.0, c="#009944",
         label="Generalized cross-validation")
plt.ylim([0.13 - 0.002, 0.17 + 0.002])
plt.yticks(range(0.13, stop=0.17, step=0.01))
plt.xlabel("Lifetime parameter \$\\lambda\$")
plt.ylabel("Loss function")
plt.legend(frameon=false)
plt.subplots_adjust(left=0.205, right=0.96, top=0.854, bottom=0.165)
plt.savefig("replication/weather/weather_gcv.png", dpi=300)



# CIs
(data, x_min, x_max, y_min, y_max) = load_data(limit=nothing)
n = nrow(data)
X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]

n_trees = 200
debias_order = 0
x_evals_original = [(20, 1020), (70, 1000), (80, 990)]
x_evals = [((x[1]-x_min)/(x_max-x_min), (x[2]-y_min)/(y_max-y_min))
           for x in x_evals_original]

forest = MondrianForest(lambda, n_trees, x_evals, 0,
                        0.05, X, Y, true, false)

for i in 1:length(x_evals)
    println()
    println("eval: ", x_evals_original[i])
    println("mu_hat: ", forest.mu_hat[i])
    println("CI: ", forest.confidence_band[i])
end
