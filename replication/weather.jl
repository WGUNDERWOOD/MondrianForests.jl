using CSV
using DataFrames
using PyPlot
using MondrianForests
using Random
using Colors
using Plots

Random.seed!(314159)
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
#rcParams["font.family"] = "Fira Sans"
#rcParams["font.sans-serif"] = "Fira Sans Compressed"
rcParams["text.usetex"] = true
plt.ioff()

function get_splits(tree::MondrianTree{d}) where {d}
    if !isnothing(tree.split_axis)
        lower = tree.tree_right.cell.lower
        upper = tree.tree_left.cell.upper
        return [(lower, upper); get_splits(tree.tree_left); get_splits(tree.tree_right)]
    else
        return Tuple{NTuple{d,Float64},NTuple{d,Float64}}[]
    end
end

function get_cells(tree::MondrianTree{d}) where {d}
    if !isnothing(tree.split_axis)
        return [get_cells(tree.tree_left); get_cells(tree.tree_right)]
    else
        return [tree.cell]
    end
end

function plot_mondrian_tree(forest::MondrianForest)
    tree = forest.trees[][]
    @assert isa(tree, MondrianTree{2})
    splits = get_splits(tree)
    cells = get_cells(tree)
    (fig, ax) = plt.subplots(figsize=(4, 3.5))

    # plot root cell
    lw = 0.3
    (l1, l2) = tree.cell.lower
    (u1, u2) = tree.cell.upper
    PyPlot.plot([l1, l1], [l1, u2], color="k", lw=lw)
    PyPlot.plot([u1, u1], [l2, u1], color="k", lw=lw)
    PyPlot.plot([l1, u1], [l2, l2], color="k", lw=lw)
    PyPlot.plot([u1, l1], [u2, u2], color="k", lw=lw)

    # define colormap
    hue1 = 25
    hue2 = 267
    colormap = ColorMap(diverging_palette(hue1, hue2, c=0.5))
    #colormap = ColorMap("RdBu")
    #colormap = cgrad([:orange, :blue], [0.1, 0.3, 0.8])
    #colormap = ColorMap(cmap("D1", reverse=true))

    # plot cells
    cell_centers = MondrianForests.get_center.(cells)
    counts = [sum(MondrianForests.are_in_same_cell(c, forest.X_data[i], tree)
                  for i in 1:forest.n_data) for c in cell_centers]
    ones = [sum(MondrianForests.are_in_same_cell(c, forest.X_data[i], tree) *
                forest.Y_data[i] for i in 1:forest.n_data) for c in cell_centers]
    colors = colormap.(ones ./ counts)
    for i in 1:length(colors)
        colors[i] = (colors[i][1], colors[i][2], colors[i][3], colors[i][4] * 0.8)
    end
    for i in 1:length(cells)
        cell = cells[i]
        x1s = [cell.lower[1], cell.lower[1], cell.upper[1], cell.upper[1]]
        x2s = [cell.lower[2], cell.upper[2], cell.upper[2], cell.lower[2]]
        fill(x1s, x2s, facecolor=colors[i])
    end

    # plot splits
    for split in splits
        x1s = [point[1] for point in split]
        x2s = [point[2] for point in split]
        PyPlot.plot(x1s, x2s, ms=0, color="k", lw=lw)
    end

    # plot data
    dry_color = "#f52e00"
    wet_color = "#005599"
    plt.scatter(data.Humidity3pm, data.Pressure3pm,
                c = [[dry_color, wet_color][Int(data.RainTomorrow[i]) + 1]
                     for i in 1:forest.n_data],
                s=2, alpha=0.3, marker=".", ec=nothing)

    # format plot
    plt.xticks([0, 1], labels=hum_lims)
    plt.yticks([0, 1], labels=pres_lims)
    plt.xlabel("Relative humidity at 3pm [%]")
    plt.ylabel("Pressure at 3pm [mbar]")
    ax.xaxis.set_label_coords(0.5, -0.04)
    ax.yaxis.set_label_coords(-0.04, 0.5)
    ax.tick_params(color="w", direction="in", pad=0)
    for side in ["bottom", "top", "left", "right"]
        ax.spines[side].set_color("#FFFFFF00")
    end

    # color key
    dry_handle = plt.scatter([], [], c=dry_color,
                             label="Dry")
    wet_handle = plt.scatter([], [], c=wet_color,
                             label="Wet")
    ax.legend([dry_handle, wet_handle], ["Dry", "Wet"],
              frameon=false,
              bbox_to_anchor=(0.95, 1))

    # color bar
    plt.colorbar(shrink=0.7, anchor=(0, 0.2), cmap=colorbar)

    return (fig, ax)
end

function normalize_data(x)
    eps = 0.001
    x_new = (x .- minimum(x)) / (maximum(x) - minimum(x))
    x_new = (1 - 2 * eps) .* x_new .+ eps
end


data = DataFrame(CSV.File("replication/weather.csv", limit=1000, missingstring="NA"))

#print(names(data))
data = data[:, [:Pressure3pm, :Humidity3pm, :MaxTemp, :RainTomorrow]]
println(nrow(data))
dropmissing!(data)
println(nrow(data))
data.RainTomorrow = replace(data.RainTomorrow, "Yes" => 1, "No" => 0)
data = Float64.(data)

hum_lims = [15, 100]
pres_lims = [1000, 1030]
data = data[hum_lims[1] .<= data.Humidity3pm .<= hum_lims[2], :]
data = data[pres_lims[1] .<= data.Pressure3pm .<= pres_lims[2], :]
display(data)
data.Humidity3pm .+= 1 * (rand(nrow(data)) .- 0.5)
data.Pressure3pm .+= 0.5 * (rand(nrow(data)) .- 0.5)
display(data)

for name in [:Pressure3pm, :Humidity3pm, :MaxTemp]
    data[:, name] = normalize_data(data[:, name])
end

lambda = 5.0
n_trees = 1
x_evals = Tuple{Float64, Float64}[]
debias_order = 0
significance_level = 0.95
estimate_var = false
#X = Array{Float64}(data[:, [:Pressure3pm, :Humidity3pm]])
#Y = Array{Float64}(data.RainTomorrow)
#Y = reshape(Y, (length(Y), 1))
X = [ntuple(j -> data[i, [:Humidity3pm, :Pressure3pm][j]], 2) for i in 1:nrow(data)]
Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]

forest = MondrianForest(lambda, n_trees, x_evals, debias_order, significance_level, X, Y, estimate_var)

#show(forest.trees[1][1])


#display(data)


#plt.close()
#plt.scatter(data.Humidity3pm, data.Pressure3pm,
            #c = data.RainTomorrow, s = 30, alpha=1, ec=nothing,
            #cmap="Paired")


#plt.scatter(data.Pressure3pm, data.MaxTemp)
#plt.scatter(data.Humidity3pm, data.MaxTemp)

(fig, ax) = plot_mondrian_tree(forest)
plt.tight_layout()
PyPlot.savefig("replication/weather.pgf")
PyPlot.savefig("replication/weather.pdf")
plt.close("all")

println()
