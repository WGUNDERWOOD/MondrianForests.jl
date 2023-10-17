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
        lower = tree.tree_right.cell.lower
        upper = tree.tree_left.cell.upper
        return [(lower, upper); get_splits(tree.tree_left); get_splits(tree.tree_right)]
    else
        return Tuple{NTuple{d,Float64},NTuple{d,Float64}}[]
    end
end

function plot_theorem_restriction(tree, cell)
    splits = get_splits(tree)
    (fig, ax) = plt.subplots(figsize=(2.2, 2.2))

    # plot root cell and extra cell
    colors = ["k", "k"]
    styles = ["solid", "dashed"]
    lw = 0.9
    for i in 1:2
        c = [tree.cell, cell][i]
        color = colors[i]
        style = styles[i]
        (l1, l2) = c.lower
        (u1, u2) = c.upper
        plt.plot([l1, l1], [l2, u2], color=color, lw=lw, linestyle=style)
        plt.plot([u1, u1], [l2, u2], color=color, lw=lw, linestyle=style)
        plt.plot([l1, u1], [l2, l2], color=color, lw=lw, linestyle=style)
        plt.plot([u1, l1], [u2, u2], color=color, lw=lw, linestyle=style)
    end

    # annotate extra cell
    plt.text(0.39, 0.33, "\$C\$", fontsize=10)

    # plot splits
    for split in splits
        x1s = [point[1] for point in split]
        x2s = [point[2] for point in split]
        plt.plot(x1s, x2s, ms=0, color="k", lw=lw)
    end

    # format plot
    plt.ylim([-0.05, 1.01])
    plt.xticks([0, 1])
    plt.yticks([0, 1])
    plt.xlabel("\$x_1\$")
    plt.ylabel("\$x_2\$")
    ax.xaxis.set_label_coords(0.5, -0.04)
    ax.yaxis.set_label_coords(-0.04, 0.5)
    ax.tick_params(color="w", direction="in", pad=0)
    ax.set_aspect("equal")
    for side in ["bottom", "top", "left", "right"]
        ax.spines[side].set_color("#FFFFFF00")
    end
    plt.tight_layout()
    return (fig, ax)
end

function plot_theorem_distribution(tree, point)
    splits = get_splits(tree)
    (fig, ax) = plt.subplots(figsize=(2.2, 2.2))

    # plot root cell
    lw = 0.9
    (l1, l2) = tree.cell.lower
    (u1, u2) = tree.cell.upper
    plt.plot([l1, l1], [l2, u2], color="k", lw=lw)
    plt.plot([u1, u1], [l2, u2], color="k", lw=lw)
    plt.plot([l1, u1], [l2, l2], color="k", lw=lw)
    plt.plot([u1, l1], [u2, u2], color="k", lw=lw)

    # get cell containing point
    cells = MondrianForests.get_cells(tree)
    cell = [c for c in cells if MondrianForests.is_in(point, c)][]

    # plot point and distribution
    e1 = 0.03
    e2 = 0.1
    plt.arrow(point[1] + e1, point[2], cell.upper[1] - point[1] - e2, 0,
              head_width=0.02, head_length=0.03, ec="k", fc="k")
    plt.arrow(point[1] - e1, point[2], cell.lower[1] - point[1] + e2, 0,
              head_width=0.02, head_length=0.03, ec="k", fc="k")
    plt.arrow(point[1], point[2] + e1, 0, cell.upper[2] - point[2] - e2,
              head_width=0.02, head_length=0.03, ec="k", fc="k")
    plt.arrow(point[1], point[2] - e1, 0, cell.lower[2] - point[2] + e2,
              head_width=0.02, head_length=0.03, ec="k", fc="k")
    plt.scatter(point[1], point[2], color="k", s=2)

    # annotate point and arrows
    plt.text(0.44, 0.28, "\$x\$", fontsize=9)
    plt.text(0.07, 0.29, "\$x_1\$", fontsize=9)
    plt.text(0.28, 0.04, "\$x_2\$", fontsize=9)
    plt.text(0.61, 0.15, "\$E_{12}\$", fontsize=9)
    plt.text(0.44, 0.65, "\$E_{22}\$", fontsize=9)

    # plot splits
    for split in splits
        x1s = [point[1] for point in split]
        x2s = [point[2] for point in split]
        plt.plot(x1s, x2s, ms=0, color="k", lw=lw)
    end

    # format plot
    plt.ylim([-0.05, 1.01])
    plt.xticks([0, 1])
    plt.yticks([0, 1])
    plt.xlabel("\$x_1\$")
    plt.ylabel("\$x_2\$")
    ax.xaxis.set_label_coords(0.5, -0.04)
    ax.yaxis.set_label_coords(-0.04, 0.5)
    ax.tick_params(color="w", direction="in", pad=0)
    ax.set_aspect("equal")
    for side in ["bottom", "top", "left", "right"]
        ax.spines[side].set_color("#FFFFFF00")
    end
    plt.tight_layout()
    return (fig, ax)
end

# construct a good tree
d = 2
lambda = 2.0
Random.seed!(0)
min_vol = 0.0
n_cells = 1
while min_vol < 0.2 || n_cells != 4
    global tree = MondrianTree(d, lambda)
    global cells = MondrianForests.get_cells(tree)
    global min_vol = minimum(MondrianForests.get_volume(c) for c in cells)
    global n_cells = length(cells)
end

# restriction theorem plot
println("plotting restriction theorem")
dpi = 500
tree = MondrianForests.restrict(tree, 1.5)
cell = MondrianCell((0.5, 0.3), (0.9, 0.85))
(fig, ax) = plot_theorem_restriction(tree, cell)
plt.savefig("replication/theorem_diagrams/theorem_restriction.png", dpi=dpi)
plt.close("all")

# distribution theorem plot
println("plotting distribution theorem")
point = (0.4, 0.25)
(fig, ax) = plot_theorem_distribution(tree, point)
plt.savefig("replication/theorem_diagrams/theorem_distribution.png", dpi=dpi)
plt.close("all")
