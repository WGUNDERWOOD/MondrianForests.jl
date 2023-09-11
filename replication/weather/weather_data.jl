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
    return (data, x_min, x_max, y_min, y_max)
end

function format_plot(ax)
    # ticks and labels
    xticks = [0, 20, 40, 60, 80, 100]
    yticks = [990 + i * 10 for i in 0:5]
    xticklabels = "\$" .* string.(xticks) .* "\$"
    yticklabels = "\$" .* string.(yticks) .* "\$"
    plt.xticks((xticks .- x_min) ./ (x_max - x_min), labels=xticklabels)
    plt.yticks((yticks .- y_min) ./ (y_max - y_min), labels=yticklabels)
    plt.xlabel("Relative humidity at 3pm (\\%)")
    plt.ylabel("Pressure at 3pm (mbar)")
    # color key
    dry_color = "#d55e00"
    wet_color = "#0035dd"
    dry_handle = plt.scatter([], [], c=dry_color)
    wet_handle = plt.scatter([], [], c=wet_color)
    ax.legend([dry_handle, wet_handle], ["Dry tomorrow", "Wet tomorrow"],
              handletextpad=0.1, frameon=false,
              bbox_to_anchor=(0.93, 1.16), ncol=2)
    # layout
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.tight_layout()
end

function plot_data()
    dry_color = "#d55e00"
    wet_color = "#0035dd"
    colors = [dry_color, wet_color][Int.(data.RainTomorrow) .+ 1]
    plt.scatter(data.Humidity3pm, data.Pressure3pm, c=colors,
                s=5, alpha=0.3, marker=".", ec=nothing)
end

function plot_splits(tree)
    # plot root cell
    splits = MondrianForests.get_splits(tree)
    lw = 0.7
    (l1, l2) = tree.cell.lower
    (u1, u2) = tree.cell.upper
    PyPlot.plot([l1, l1], [l1, u2], color="k", lw=lw)
    PyPlot.plot([u1, u1], [l2, u1], color="k", lw=lw)
    PyPlot.plot([l1, u1], [l2, l2], color="k", lw=lw)
    PyPlot.plot([u1, l1], [u2, u2], color="k", lw=lw)
    # plot splits
    for split in splits
        x1s = [point[1] for point in split]
        x2s = [point[2] for point in split]
        PyPlot.plot(x1s, x2s, ms=0, color="k", lw=lw)
    end
end

function plot_cells(tree)
    cells = MondrianForests.get_cells(tree)
    X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
    Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]
    n = length(X)
    cell_centers = MondrianForests.get_center.(cells)
    counts = [sum(MondrianForests.are_in_same_cell(c, X[i], tree)
                  for i in 1:n) for c in cell_centers]
    ones = [sum(MondrianForests.are_in_same_cell(c, X[i], tree) *
                Y[i] for i in 1:n) for c in cell_centers]
    hue1 = 25
    hue2 = 267
    colormap = ColorMap(diverging_palette(hue1, hue2, c=0.5))
    colors = colormap.(ones ./ counts)
    for i in 1:length(cells)
        cell = cells[i]
        x1s = [cell.lower[1], cell.lower[1], cell.upper[1], cell.upper[1]]
        x2s = [cell.lower[2], cell.upper[2], cell.upper[2], cell.lower[2]]
        fill(x1s, x2s, facecolor=colors[i])
    end
end

function plot_small_forest(trees, ax)

    # get all cells
    all_cells = [MondrianForests.get_cells(tree) for tree in trees]

    # get refined cells
    refined_cells = copy(all_cells[1])
    if length(trees) >= 2
        for j in 2:length(trees)
            for cell in all_cells[j]
                new_cells = [MondrianForests.get_intersection(cell, c) for c in refined_cells]
                new_cells = [c for c in new_cells if !isnothing(c)]
                append!(refined_cells, new_cells)
                refined_cells = unique(refined_cells)
            end
        end
    end

    # get counts
    X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
    Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]
    n = length(X)

    ratios = Float64[]

    for r_cell in refined_cells
        x = MondrianForests.get_center(r_cell)
        cells = [[cell for cell in all_cells[j] if
                  MondrianForests.is_in(x, cell)][] for j in 1:length(trees)]
        counts = [sum(MondrianForests.is_in(X[i], cell) for i in 1:n) for cell in cells]
        ones = [sum(MondrianForests.is_in(X[i], cell) * Y[i] for i in 1:n) for cell in cells]
        ratio = sum(ones ./ counts) / length(trees)
        push!(ratios, ratio)
    end
    hue1 = 25
    hue2 = 267
    colormap = ColorMap(diverging_palette(hue1, hue2, c=0.5))
    colors = colormap.(ratios)
    for i in 1:length(refined_cells)
        cell = refined_cells[i]
        x1s = [cell.lower[1], cell.lower[1], cell.upper[1], cell.upper[1]]
        x2s = [cell.lower[2], cell.upper[2], cell.upper[2], cell.lower[2]]
        fill(x1s, x2s, facecolor=colors[i])
    end
end

function plot_large_forest(forest, ax)
    n_x1s = length(x1s)
    n_x2s = length(x2s)
    reshaped_mu_hat = reshape(forest.mu_hat, (n_x2s, n_x1s))
    hue1 = 25
    hue2 = 267
    colormap = ColorMap(diverging_palette(hue1, hue2, c=0.5))
    plt.contourf(x1s, x2s, reshaped_mu_hat, cmap=colormap, levels=1000)

end

function make_data_plot(data, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=(4, 3.5))
    plot_data()
    format_plot(ax)
    PyPlot.savefig(filename, dpi=300)
    plt.close("all")
end

function make_data_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=(4, 3.5))
    plot_splits(tree)
    plot_data()
    format_plot(ax)
    PyPlot.savefig(filename, dpi=300)
    plt.close("all")
end

function make_data_filled_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=(4, 3.5))
    plot_splits(tree)
    plot_cells(tree)
    plot_data()
    format_plot(ax)
    PyPlot.savefig(filename, dpi=300)
    plt.close("all")
end

function make_filled_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=(4, 3.5))
    plot_splits(tree)
    plot_cells(tree)
    format_plot(ax)
    PyPlot.savefig(filename, dpi=300)
    plt.close("all")
end

function make_small_forest_plot(data, trees, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=(4, 3.5))
    plot_small_forest(trees, ax)
    format_plot(ax)
    PyPlot.savefig(filename, dpi=300)
    plt.close("all")
end

function make_large_forest_plot(data, forest, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=(4, 3.5))
    plot_large_forest(forest, ax)
    format_plot(ax)
    PyPlot.savefig(filename, dpi=300)
    plt.close("all")
end

(data, x_min, x_max, y_min, y_max) = load_data(limit=nothing)
lambda = 5.0

# plot data
println("plotting data")
filename = "replication/weather_data.png"
#make_data_plot(data, x_min, x_max, y_min, y_max, filename)

# plot data and partition
println("plotting data and partition")
Random.seed!(314159)
tree = MondrianTree(2, lambda)
filename = "replication/weather_data_partition.png"
#make_data_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)

# plot data and filled partition
println("plotting data and filled partition")
filename = "replication/weather_data_filled_partition.png"
#make_data_filled_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)

# plot filled partition
println("plotting filled partition")
filename = "replication/weather_filled_partition.png"
#make_filled_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)

seeds = [314160, 314161, 314162, 314163, 314164, 314165]
trees = []
for i in 1:length(seeds)
    println("plotting filled partition for tree " * string(i))
    global seed = seeds[i]
    Random.seed!(seed)
    tree = MondrianTree(2, lambda)
    push!(trees, tree)
    global filename = "replication/weather_filled_partition_" * string(i) * ".png"
    #make_filled_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)
end

# plot small forest
debias_order = 0
significance_level = 0.95
estimate_var = false
n_trees = 2
X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]
x_evals = Tuple{Float64, Float64}[]
forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                        significance_level, X, Y, estimate_var)

filename = "replication/weather_forest_small.png"
make_small_forest_plot(data, trees, x_min, x_max, y_min, y_max, filename)

# plot large forest
#n_x1_evals = 50
#n_x2_evals = 50
#x1s = range(1 / n_x1_evals, 1 - 1 / n_x1_evals, length=n_x1_evals)
#x2s = range(1 / n_x2_evals, 1 - 1 / n_x2_evals, length=n_x2_evals)
#x_evals = Tuple{Float64, Float64}[(x1s[i], x2s[j]) for i in 1:n_x1_evals for j in 1:n_x2_evals]

#n_trees_all = [1, 2, 3, 10, 100]
#for i in 1:length(n_trees_all)
    #n_trees = n_trees_all[i]
    #println("computing forest with ", n_trees, " trees")
    #seed = 314165
    #Random.seed!(seed)
    #forest = MondrianForest(lambda, n_trees, x_evals, debias_order,
                            #significance_level, X, Y, estimate_var)

    #println("plotting forest with ", n_trees, " trees")
    #global filename = "replication/weather_forest_" * string(i) * ".png"
    #make_forest_plot(data, forest, x_min, x_max, y_min, y_max, filename)
#end

