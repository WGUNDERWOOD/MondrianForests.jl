using CSV
using DataFrames
using PyPlot
using Random
using Colors
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
    Random.seed!(1)
    data = shuffle(data)
    println("number of samples ", nrow(data))
    return (data, x_min, x_max, y_min, y_max)
end

function format_plot(ax)
    # ticks and labels
    xticks = [0, 20, 40, 60, 80, 100]
    yticks = [990 + i * 10 for i in 0:5]
    xticklabels = "\$" .* string.(xticks) .* "\$"
    yticklabels = "\$" .* string.(yticks) .* "\$"
    plt.xticks((xticks .- x_min) ./ (x_max - x_min), labels=xticklabels, fontsize=11)
    plt.yticks((yticks .- y_min) ./ (y_max - y_min), labels=yticklabels, fontsize=11)
    plt.xlabel("Relative humidity at 3pm (\\%)", fontsize=12)
    plt.ylabel("Pressure at 3pm (mbar)", fontsize=12)
    # color key
    #handle = plt.scatter([], [], c="white")
    dry_handle = plt.scatter([], [], c=dry_color)
    wet_handle = plt.scatter([], [], c=wet_color)
    ax.legend([wet_handle, dry_handle],
              ["Rain next day", "No rain next day"],
              handletextpad=0.1, frameon=false,
              bbox_to_anchor=(1.04, 1.13), ncol=3,
              fontsize=12,
              columnspacing=0.8)
    # layout
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    return plt.tight_layout()
end

function plot_data()
    colors = [dry_color, wet_color][Int.(data.RainTomorrow) .+ 1]
    return plt.scatter(data.Humidity3pm, data.Pressure3pm, c=colors,
                       s=1, alpha=0.3, marker=".", ec=nothing)
end

function plot_splits(tree)
    # plot root cell
    split_trees = [t for t in get_subtrees(tree) if t.is_split]
    splits = [(t.tree_right.lower, t.tree_left.upper) for t in split_trees]
    lw = 0.7
    (l1, l2) = tree.lower
    (u1, u2) = tree.upper
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

function get_cell_colormap()
    return ColorMap(diverging_palette(27, 225; c=0.6, s=0.8))
end

function plot_cells(tree)
    cells = MondrianForests.get_leaves(tree)
    X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
    Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]
    n = length(X)
    cell_centers = MondrianForests.get_center.(cells)
    counts = [sum(MondrianForests.are_in_same_leaf(c, X[i], tree)
                  for i in 1:n) for c in cell_centers]
    ones = [sum(MondrianForests.are_in_same_leaf(c, X[i], tree) *
                Y[i] for i in 1:n) for c in cell_centers]
    colormap = get_cell_colormap()
    colors = colormap.(ones ./ counts)
    for i in 1:length(cells)
        cell = cells[i]
        x1s = [cell.lower[1], cell.lower[1], cell.upper[1], cell.upper[1]]
        x2s = [cell.lower[2], cell.upper[2], cell.upper[2], cell.lower[2]]
        fill(x1s, x2s, facecolor=colors[i])
    end
end

function plot_forest(trees, ax)
    # get all cells and refinement
    refined_tree = get_common_refinement(trees)
    refined_cells = get_leaves(refined_tree)
    all_cells = [get_leaves(t) for t in trees]
    # get counts
    X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
    Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]
    n = length(X)
    all_counts = [[sum(MondrianForests.is_in(X[i], cell) for i in 1:n)
                   for cell in all_cells[j]] for j in 1:length(trees)]
    all_ones = [[sum(MondrianForests.is_in(X[i], cell) * Y[i] for i in 1:n)
                 for cell in all_cells[j]] for j in 1:length(trees)]
    ratios = zeros(length(refined_cells))

    for c in 1:length(refined_cells)
        r_cell = refined_cells[c]
        x = MondrianForests.get_center(r_cell)
        cell_ids = [[c for c in 1:length(all_cells[j]) if
                     MondrianForests.is_in(x, all_cells[j][c])][] for j in 1:length(trees)]
        counts = [all_counts[j][cell_ids[j]] for j in 1:length(trees)]
        ones = [all_ones[j][cell_ids[j]] for j in 1:length(trees)]
        rs = ones ./ counts
        if !all(isnan.(rs))
            ratios[c] = sum(r for r in rs if !isnan(r)) / sum(1 for r in rs if !isnan(r))
        else
            ratios[c] = NaN
        end
    end

    # plot
    colormap = get_cell_colormap()
    colors = colormap.(ratios)
    for i in 1:length(refined_cells)
        println(i, " / ", length(refined_cells))
        cell = refined_cells[i]
        x1s = [cell.lower[1], cell.lower[1], cell.upper[1], cell.upper[1]]
        x2s = [cell.lower[2], cell.upper[2], cell.upper[2], cell.lower[2]]
        fill(x1s, x2s, facecolor=colors[i])
    end
end

function plot_debiased_forest(trees1, trees2, ax)
    # get all cells and refinement
    refined_tree = get_common_refinement([trees1; trees2])
    refined_cells = get_leaves(refined_tree)
    all_cells1 = [get_leaves(t) for t in trees1]
    all_cells2 = [get_leaves(t) for t in trees2]
    # get counts
    X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
    Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]
    n = length(X)
    all_counts1 = [[sum(MondrianForests.is_in(X[i], cell) for i in 1:n)
                    for cell in all_cells1[j]] for j in 1:length(trees1)]
    all_counts2 = [[sum(MondrianForests.is_in(X[i], cell) for i in 1:n)
                    for cell in all_cells2[j]] for j in 1:length(trees2)]
    all_ones1 = [[sum(MondrianForests.is_in(X[i], cell) * Y[i] for i in 1:n)
                  for cell in all_cells1[j]] for j in 1:length(trees1)]
    all_ones2 = [[sum(MondrianForests.is_in(X[i], cell) * Y[i] for i in 1:n)
                  for cell in all_cells2[j]] for j in 1:length(trees2)]
    ratios = zeros(length(refined_cells))
    ratios1 = zeros(length(refined_cells))
    ratios2 = zeros(length(refined_cells))
    debias_order = 1
    debias_coeffs = MondrianForests.get_debias_coeffs(debias_order)

    for c in 1:length(refined_cells)
        r_cell = refined_cells[c]
        x = MondrianForests.get_center(r_cell)
        cell_ids1 = [[c
                      for c in 1:length(all_cells1[j])
                      if
                      MondrianForests.is_in(x, all_cells1[j][c])][] for j in 1:length(trees1)]
        cell_ids2 = [[c
                      for c in 1:length(all_cells2[j])
                      if
                      MondrianForests.is_in(x, all_cells2[j][c])][] for j in 1:length(trees2)]
        counts1 = [all_counts1[j][cell_ids1[j]] for j in 1:length(trees1)]
        counts2 = [all_counts2[j][cell_ids2[j]] for j in 1:length(trees2)]
        ones1 = [all_ones1[j][cell_ids1[j]] for j in 1:length(trees1)]
        ones2 = [all_ones2[j][cell_ids2[j]] for j in 1:length(trees2)]
        rs1 = ones1 ./ counts1
        rs2 = ones2 ./ counts2
        if !all(isnan.(rs1))
            ratios1[c] = sum(r for r in rs1 if !isnan(r)) / sum(1 for r in rs1 if !isnan(r))
        else
            ratios1[c] = NaN
        end
        if !all(isnan.(rs2))
            ratios2[c] = sum(r for r in rs2 if !isnan(r)) / sum(1 for r in rs2 if !isnan(r))
        else
            ratios2[c] = NaN
        end

        if !isnan(ratios1[c]) && !isnan(ratios2[c])
            ratios[c] = ratios1[c] * debias_coeffs[1] + ratios2[c] * debias_coeffs[2]
        elseif !isnan(ratios1[c])
            ratios[c] = ratios1[c]
        elseif !isnan(ratios2[c])
            ratios[c] = ratios2[c]
        else
            ratios[c] = NaN
        end
    end

    # plot
    colormap = get_cell_colormap()
    colors = colormap.(ratios)
    for i in 1:length(refined_cells)
        println(i, " / ", length(refined_cells))
        cell = refined_cells[i]
        x1s = [cell.lower[1], cell.lower[1], cell.upper[1], cell.upper[1]]
        x2s = [cell.lower[2], cell.upper[2], cell.upper[2], cell.lower[2]]
        fill(x1s, x2s, facecolor=colors[i])
    end
end

function make_data_plot(data, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=figsize)
    plot_data()
    format_plot(ax)
    PyPlot.savefig(filename, dpi=dpi)
    return plt.close("all")
end

function make_data_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=figsize)
    plot_splits(tree)
    plot_data()
    format_plot(ax)
    PyPlot.savefig(filename, dpi=dpi)
    return plt.close("all")
end

function make_data_filled_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=figsize)
    plot_splits(tree)
    plot_cells(tree)
    plot_data()
    format_plot(ax)
    PyPlot.savefig(filename, dpi=dpi)
    return plt.close("all")
end

function make_filled_partition_plot(data, tree, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=figsize)
    plot_splits(tree)
    plot_cells(tree)
    format_plot(ax)
    PyPlot.savefig(filename, dpi=dpi)
    return plt.close("all")
end

function make_forest_plot(data, trees, x_min, x_max, y_min, y_max, filename)
    (fig, ax) = plt.subplots(figsize=figsize)
    plot_forest(trees, ax)
    format_plot(ax)
    PyPlot.savefig(filename, dpi=dpi)
    return plt.close("all")
end

function make_forest_design_plot(data, trees, x_min, x_max,
                                 y_min, y_max, design_points, filename)
    (fig, ax) = plt.subplots(figsize=figsize)
    plot_forest(trees, ax)
    for i in 1:length(design_points)
        p = design_points[i]
        plt.scatter(p[1], p[2], c="k", s=30, marker="x", linewidths=1)
        plt.text(p[1] + 0.07, p[2] - 0.003, "\$$(string(i))\$", fontsize=12,
                 ha="center", va="center")
    end
    format_plot(ax)
    PyPlot.savefig(filename, dpi=dpi)
    return plt.close("all")
end

function make_debiased_forest_design_plot(data, trees1, trees2, x_min, x_max,
                                          y_min, y_max, design_points, filename)
    (fig, ax) = plt.subplots(figsize=figsize)
    plot_debiased_forest(trees1, trees2, ax)
    for i in 1:length(design_points)
        p = design_points[i]
        plt.scatter(p[1], p[2], c="k", s=30, marker="x", linewidths=1)
        plt.text(p[1] + 0.07, p[2] - 0.003, "\$$(string(i))\$", fontsize=12,
                 ha="center", va="center")
    end
    format_plot(ax)
    PyPlot.savefig(filename, dpi=dpi)
    return plt.close("all")
end

# get data and plot params
limit = nothing
(data, x_min, x_max, y_min, y_max) = load_data(limit=limit)
dry_color = "#da6200"
wet_color = "#0080d0"
figsize = (4, 4.2)
dpi = 350

# make trees
seeds = 4:54
trees = MondrianTree{2}[]
lambda = 5.0
for i in 1:length(seeds)
    seed = seeds[i]
    Random.seed!(seed)

    if i <= 5
        global min_vol = 0.0
        while min_vol < 0.009
            tree = MondrianTree(2, lambda)
            global cells = MondrianForests.get_leaves(tree)
            global min_vol = minimum(MondrianForests.get_volume(c) for c in cells)
        end
    else
        tree = MondrianTree(2, lambda)
    end

    push!(trees, tree)
end

# plot data
println("plotting data")
filename = "./replication/weather/weather_data.png"
make_data_plot(data, x_min, x_max, y_min, y_max, filename)

# plot data and filled partition
println("plotting data and filled partition")
filename = "./replication/weather/weather_data_filled_partition.png"
make_data_filled_partition_plot(data, trees[1], x_min, x_max, y_min, y_max, filename)

# plot small forest
significance_level = 0.95
estimate_var = false
n_trees = 2
X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]
x_evals = Tuple{Float64,Float64}[]
println("plotting forest with ", n_trees, " trees")
global filename = "./replication/weather/weather_forest_" * string(n_trees) * ".png"
make_forest_plot(data, trees[1:n_trees], x_min, x_max, y_min, y_max, filename)

# plot debiased forest with design points
n_trees = 20
debias_order = 1
debias_scaling = MondrianForests.get_debias_scaling(debias_order)
lambda2 = lambda * debias_scaling[2]
trees2 = [MondrianTree(2, lambda2) for _ in 1:n_trees]
println("plotting debiased forest with ", n_trees, " trees and design points")
design_points = [(20, 1020), (70, 1000), (80, 990)]
design_points = [((x[1] - x_min) / (x_max - x_min), (x[2] - y_min) / (y_max - y_min))
                 for x in design_points]
global filename = "./replication/weather/weather_debiased_forest_design.png"
make_debiased_forest_design_plot(data, trees[1:n_trees], trees2[1:n_trees],
                                 x_min, x_max, y_min, y_max, design_points,
                                 filename)

# plot forest with design points
n_trees = 40
println("plotting forest with ", n_trees, " trees and design points")
global filename = "./replication/weather/weather_forest_design.png"
make_forest_design_plot(data, trees[1:n_trees], x_min, x_max, y_min,
                        y_max, design_points, filename)
