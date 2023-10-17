using MondrianForests
using PyPlot
using Random

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

function plot_mondrian_tree(tree::MondrianTree)
    @assert isa(tree, MondrianTree{2})
    splits = get_splits(tree)
    (fig, ax) = plt.subplots(figsize=(2.1, 2.1))

    # plot root cell
    lw = 0.3
    (l1, l2) = tree.cell.lower
    (u1, u2) = tree.cell.upper
    plot([l1, l1], [l1, u2], color="k", lw=lw)
    plot([u1, u1], [l2, u1], color="k", lw=lw)
    plot([l1, u1], [l2, l2], color="k", lw=lw)
    plot([u1, l1], [u2, u2], color="k", lw=lw)

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

d = 2
lambdas = [3.0, 10.0, 30.0]
Random.seed!(314159)
git_root = strip(read(`git rev-parse --show-toplevel`, String), '\n')

for i in 1:length(lambdas)
    lambda = lambdas[i]
    tree = MondrianTree(d, lambda)
    (fig, ax) = plot_mondrian_tree(tree)
    savefig("./replication/partition_plots/plot_mondrian_process_$i.pgf", bbox_inches="tight")
    savefig("./replication/partition_plots/plot_mondrian_process_$i.pdf", bbox_inches="tight")
    plt.close("all")
end
