using MondrianForests
using PyPlot
using Random
using Revise

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

function plot_mondrian_tree(tree::MondrianTree, cell::Union{MondrianCell, Nothing})

    splits = get_splits(tree)
    (fig, ax) = plt.subplots(figsize=(2.1, 2.1))

    # plot root cell
    lw = 0.9
    (l1, l2) = tree.cell.lower
    (u1, u2) = tree.cell.upper
    plot([l1, l1], [l2, u2], color="k", lw=lw)
    plot([u1, u1], [l2, u2], color="k", lw=lw)
    plot([l1, u1], [l2, l2], color="k", lw=lw)
    plot([u1, l1], [u2, u2], color="k", lw=lw)

    # plot splits
    for split in splits
        x1s = [point[1] for point in split]
        x2s = [point[2] for point in split]
        plot(x1s, x2s, ms=0, color="k", lw=lw)
    end

    # highlight cell
    if !isnothing(cell)
        (l1, l2) = cell.lower
        (u1, u2) = cell.upper
        plot([l1, l1], [l2, u2], color="r", lw=lw, zorder=10)
        plot([l1, u1], [l2, l2], color="r", lw=lw, zorder=10)
    end

    # annotate cells
    ids = MondrianForests.get_ids(tree)
    cells = MondrianForests.get_cells(tree)
    centers = MondrianForests.get_center.(cells)
    for i in 1:length(cells)
        label = "\$C_{\\mathrm{$(ids[i])}}\$"
        plt.text(centers[i][1], centers[i][2], label, ha="center",
                 va="center", fontsize=10)
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

d = 2
lambda = 2.0

Random.seed!(0)
min_vol = 0.0
n_cells = 1
while min_vol < 0.2 || n_cells != 4
    global tree = MondrianTree(d, lambda)
    global cells = MondrianForests.get_cells(tree)
    global min_vol =  minimum(MondrianForests.get_volume(c) for c in cells)
    global n_cells =  length(cells)
end


println(min_vol)
println(n_cells)

times = sort(unique(MondrianForests.get_split_times(tree)))
show(tree)

for i in 1:length(times)
    println(i)
    t = times[i]
    restricted_tree = MondrianForests.restrict(tree, t)
    restricted_cells = MondrianForests.get_cells(restricted_tree)
    (fig, ax) = plot_mondrian_tree(restricted_tree, nothing)
    savefig("replication/mondrian_process/mondrian_process_$(i).png",
            bbox_inches="tight", dpi=200)
    for j in 1:length(restricted_cells)
        cell = restricted_cells[j]
        (fig, ax) = plot_mondrian_tree(restricted_tree, cell)
        savefig("replication/mondrian_process/mondrian_process_$(i)_$(j).png",
                bbox_inches="tight", dpi=200)
    end
    plt.close("all")
end
