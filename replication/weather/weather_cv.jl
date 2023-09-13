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

n_evals = 100
lambdas = range(0.5, stop=10.0, step=0.5)
n_trees = 50
debias_order = 0

eval_ids = sort(shuffle(1:n)[1:n_evals])
x_evals = X[eval_ids]
y_evals = Y[eval_ids]
#forest = MondrianForest(2.0, 10, x_evals, 0,
                        #0.05, X, Y, false, true)

#forest.gcv_dof

#=
Random.seed!()
mses = Float64[]
for lambda in lambdas
    println("lambda: ", lambda)
    mse = get_mse(n_evals, X, Y, lambda, n_trees, debias_order, n_folds)
    push!(mses, mse)
end
=#

Random.seed!()
gcvs = Float64[]
mses = Float64[]
for lambda in lambdas
    println("lambda: ", lambda)
    forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                            0.05, X, Y, false, true)
    gcv_dof = forest.gcv_dof
    println("gcv dof: ", gcv_dof)
    mse = sum((y_evals .- forest.mu_hat).^2) / n_evals
    gcv = mse / ((1 - gcv_dof / n)^2)
    push!(gcvs, gcv)
    push!(mses, mse)
end

gcvs

(fig, ax) = plt.subplots(figsize=(5, 3))
plt.plot(lambdas, gcvs)
plt.plot(lambdas, mses)
plt.savefig("replication/weather/weather_gcv.png", dpi=300)
