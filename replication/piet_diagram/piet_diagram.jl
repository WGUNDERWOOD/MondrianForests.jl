using MondrianForests
using PyPlot
using Random
using Revise

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["text.latex.preamble"] = "\\usepackage[sfdefault,light]{FiraSans}"
plt.ioff()

function get_splits(tree::MondrianTree{d}) where {d}
    if !isnothing(tree.split_axis)
        lower = tree.tree_right.lower
        upper = tree.tree_left.upper
        return [(lower, upper); get_splits(tree.tree_left); get_splits(tree.tree_right)]
    else
        return Tuple{NTuple{d,Float64},NTuple{d,Float64}}[]
    end
end

function plot_piet_mondrian(tree)
    splits = get_splits(tree)
    (fig, ax) = plt.subplots(figsize=(2.2, 2.2))

    # plot root cell
    lw = 0.9
    (l1, l2) = tree.lower
    (u1, u2) = tree.upper
    plt.plot([l1, l1], [l2, u2], color="k", lw=lw)
    plt.plot([u1, u1], [l2, u2], color="k", lw=lw)
    plt.plot([l1, u1], [l2, l2], color="k", lw=lw)
    plt.plot([u1, l1], [u2, u2], color="k", lw=lw)

    # plot splits
    for split in splits
        x1s = [point[1] for point in split]
        x2s = [point[2] for point in split]
        plt.plot(x1s, x2s, ms=0, color="k", lw=lw)
    end

    # format plot
    plt.ylim([-0.01, 1.01])
    plt.xticks([])
    plt.yticks([])
    ax.tick_params(color="w", direction="in", pad=0)
    ax.set_aspect("equal")
    for side in ["bottom", "top", "left", "right"]
        ax.spines[side].set_color("#FFFFFF00")
    end
    plt.tight_layout()
    return (fig, ax)
end

# piet plot
Random.seed!(0)
println("plotting piet")
lambda = 4.0
min_vol = 0.0
best_min_vol = 0.0
dpi = 500
while min_vol < 0.03 || n_cells < 9
    global tree = MondrianTree(2, lambda)
    global cells = get_leaves(tree)
    global min_vol = minimum(MondrianForests.get_volume(c) for c in cells)
    global n_cells = length(cells)
    if min_vol > best_min_vol
        global best_min_vol = min_vol
    end
end
(fig, ax) = plot_piet_mondrian(tree)
plt.savefig("replication/piet_diagram/piet_diagram.png", dpi=dpi,
            bbox_inches="tight", pad_inches=-0.008)
plt.close("all")
