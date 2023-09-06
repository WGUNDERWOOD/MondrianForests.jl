using CSV
using DataFrames
using PyPlot
using MondrianForests
using Random

#Random.seed!(314159)
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
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
        return [get_splits(tree.tree_left); get_splits(tree.tree_right)]
    else
        return [(tree.cell.lower, tree.cell.upper)]
    end
end

function plot_mondrian_tree(forest::MondrianForest)
    tree = forest.trees[][]
    @assert isa(tree, MondrianTree{2})
    splits = get_splits(tree)
    cells = get_cells(tree)
    (fig, ax) = plt.subplots(figsize=(2.1, 2.1))

    # plot root cell
    lw = 0.3
    (l1, l2) = tree.cell.lower
    (u1, u2) = tree.cell.upper
    plot([l1, l1], [l1, u2], color="k", lw=lw)
    plot([u1, u1], [l2, u1], color="k", lw=lw)
    plot([l1, u1], [l2, l2], color="k", lw=lw)
    plot([u1, l1], [u2, u2], color="k", lw=lw)

    # plot cells
    # TODO
    for cell in cells
        x1s = [point[1] for point in cell]
        x2s = [point[2] for point in cell]
        fill(x1s, x2s, color="k")
    end

    # plot splits
    for split in splits
        x1s = [point[1] for point in split]
        x2s = [point[2] for point in split]
        plot(x1s, x2s, ms=0, color="k", lw=lw)
    end

    # format plot
    plt.xticks([0, 1])
    plt.yticks([0, 1])
    plt.xlabel("\$x_1\$")
    plt.ylabel("\$x_2\$")
    ax.xaxis.set_label_coords(0.5, -0.04)
    ax.yaxis.set_label_coords(-0.04, 0.5)
    ax.tick_params(color="w", direction="in", pad=0)
    for side in ["bottom", "top", "left", "right"]
        ax.spines[side].set_color("#FFFFFF00")
    end

    return (fig, ax)
end


data = DataFrame(CSV.File("replication/weather.csv", limit=1000, missingstring="NA"))

#print(names(data))
data = data[:, [:Pressure3pm, :Humidity3pm, :MaxTemp, :RainTomorrow]]
println(nrow(data))
dropmissing!(data)
println(nrow(data))
data.RainTomorrow = replace(data.RainTomorrow, "Yes" => 1, "No" => 0)
data = Float64.(data)

lambda = 5.0
n_trees = 1
x_evals = Tuple{Float64, Float64}[]
debias_order = 0
significance_level = 0.95
estimate_var = false
#X = Array{Float64}(data[:, [:Pressure3pm, :Humidity3pm]])
#Y = Array{Float64}(data.RainTomorrow)
#Y = reshape(Y, (length(Y), 1))
X = [ntuple(j -> data[i, [:Pressure3pm, :Humidity3pm][j]], 2) for i in 1:nrow(data)]
Y = [data[i, :RainTomorrow] for i in 1:nrow(data)]

forest = MondrianForest(lambda, n_trees, x_evals, debias_order, significance_level, X, Y, estimate_var)

show(forest.trees[1][1])


#data = data[15 .<= data.Humidity3pm, :]
#data = data[1000 .<= data.Pressure3pm .<= 1030, :]
#display(data)

#plt.close()
#plt.scatter(data.Humidity3pm, data.Pressure3pm,
            #c = data.RainTomorrow, s = 30, alpha=1, ec=nothing,
            #cmap="Paired")


#plt.scatter(data.Pressure3pm, data.MaxTemp)
#plt.scatter(data.Humidity3pm, data.MaxTemp)

(fig, ax) = plot_mondrian_tree(forest)
savefig("weather_tree.pdf", bbox_inches="tight")
plt.close("all")

println()
