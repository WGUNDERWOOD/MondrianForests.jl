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

function plot_mondrian_process(partition)
    tree = partition["tree"]
    splits = get_splits(tree)
    (fig, ax) = plt.subplots(figsize=(2.2, 2.2))

    # highlight current cell
    if !isnothing(partition["current"])
        cell = partition["current"]
        x1s = [cell.lower[1], cell.lower[1], cell.upper[1], cell.upper[1]]
        x2s = [cell.lower[2], cell.upper[2], cell.upper[2], cell.lower[2]]
        fill(x1s, x2s, facecolor="#ecd9ff")
    end

    # highlight leaves
    for cell in partition["terminals"]
        x1s = [cell.lower[1], cell.lower[1], cell.upper[1], cell.upper[1]]
        x2s = [cell.lower[2], cell.upper[2], cell.upper[2], cell.lower[2]]
        fill(x1s, x2s, facecolor="#b5fdc7")
    end

    # plot root cell
    lw = 0.9
    (l1, l2) = tree.cell.lower
    (u1, u2) = tree.cell.upper
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

    # annotate cells
    ids = MondrianForests.get_ids(tree)
    cells = MondrianForests.get_cells(tree)
    centers = MondrianForests.get_center.(cells)
    for i in 1:length(cells)
        if ids[i] == ""
            label = "\$C_{\\emptyset}\$"
        else
            label = "\$C_{\\mathrm{$(ids[i])}}\$"
        end
        plt.text(centers[i][1], centers[i][2], label, ha="center",
                 va="center", fontsize=10)
    end

    # add split point
    subtrees = MondrianForests.get_subtrees(tree)
    current = partition["current"]
    if !isnothing(current)
        if !(current in cells)
            subtree = [c for c in subtrees if c.cell == current][]
            J = subtree.split_axis
            S = subtree.split_location
            J == 1 ? x = S - 0.005 : x = -0.06
            J == 1 ? y = -0.07 : y = S - 0.005
            plt.text(x, y, "\$S\$", fontsize=10, ha="center", va="center")
        end
    end

    # format plot
    plt.ylim([-0.05, 1.01])
    plt.xticks([0, 1])
    plt.yticks([0, 1])
    plt.xlabel("\$x_1\$")
    plt.ylabel("\$x_2\$")
    ax.xaxis.set_label_coords(0.5, -0.04)
    ax.yaxis.set_label_coords(-0.06, 0.5)
    ax.tick_params(color="w", direction="in", pad=0)
    ax.set_aspect("equal")
    for side in ["bottom", "top", "left", "right"]
        ax.spines[side].set_color("#FFFFFF00")
    end

    # color key
    current_handle = plt.scatter([], [], c="#d9b3ff")
    terminal_handle = plt.scatter([], [], c="#6feb8e")
    if length(cells) >= 4 || length(partition["terminals"]) >= 1
        ax.legend([current_handle, terminal_handle],
                  ["Current", "Leaf"], ncol=2,
                  handletextpad=0.1, frameon=false, columnspacing=0.8,
                  bbox_to_anchor=(0.47, 1.19), loc="upper center")
    elseif length(cells) >= 2 || !isnothing(partition["current"])
        ax.legend([current_handle], ["Current"], ncol=2,
                  handletextpad=0.1, frameon=false, columnspacing=0.8,
                  bbox_to_anchor=(0.277, 1.19), loc="upper center")
    end

    return (fig, ax)
end

function get_tree_info(tree::MondrianTree)
    info = (tree.id, tree.creation_time, tree.cell, tree.id, tree.split_axis)
    if !isnothing(tree.split_axis)
        return [info; get_tree_info(tree.tree_left); get_tree_info(tree.tree_right)]
    else
        return [info]
    end
end

function get_horizontal_value(id::String)
    value = 0.0
    for i in 1:length(id)
        c = id[i]
        if c == 'L'
            value -= 0.5^i
        elseif c == 'R'
            value += 0.5^i
        end
    end
    return value
end

function plot_mondrian_tree(partition)
    tree = partition["tree"]
    (fig, ax) = plt.subplots(figsize=(2.2, 2.2))
    info = get_tree_info(tree)
    ids = [i[1] for i in info]
    times = [i[2] for i in info]
    cells = [i[3] for i in info]
    n = length(times)

    # plot split points
    for i in 1:n
        if cells[i] == partition["current"]
            color = "#ecd9ff"
        elseif cells[i] in partition["terminals"]
            color = "#b5fdc7"
        else
            color = "white"
        end
        circle1 = plt.Circle((x_locs[ids[i]], times[i]), 0.2, color="k", zorder=10)
        circle2 = plt.Circle((x_locs[ids[i]], times[i]), 0.19, color=color, zorder=20)
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        if ids[i] == ""
            label = "\$\\emptyset\$"
        else
            label = "\$\\mathrm{$(ids[i])}\$"
        end
        plt.text(x_locs[ids[i]] + 0.005, times[i] + 0.01, label, ha="center",
                 va="center", fontsize=8, zorder=30)
    end

    # plot tree
    lw = 0.9
    for i in 1:n
        id = ids[i]
        if id * "L" in ids
            x_left = [x_locs[ids[i]] for i in 1:n if ids[i] == id * "L"][]
            x_right = [x_locs[ids[i]] for i in 1:n if ids[i] == id * "R"][]
            time1 = times[i]
            time2 = [times[i] for i in 1:n if ids[i] == id * "L"][]
            plt.plot([x_left, x_right], [time1, time1], color="k", lw=lw)
            plt.plot([x_left, x_left], [time1, time2], color="k", lw=lw)
            plt.plot([x_right, x_right], [time1, time2], color="k", lw=lw)
        end
    end

    # add split time
    cells = MondrianForests.get_cells(tree)
    subtrees = MondrianForests.get_subtrees(tree)
    current = partition["current"]
    if !isnothing(current)
        if !(current in cells)
            subtree = [c for c in subtrees if c.cell == current][]
            t = subtree.tree_left.creation_time
            x_left = x_locs[subtree.tree_left.id]
            x_right = x_locs[subtree.tree_right.id]
            if x_left > 1
                plt.text(x_left - 0.8, t, "\$t + E\$", fontsize=10,
                         ha="center", va="center")
                plt.plot([x_left - 0.47, 3], [t, t], linestyle="dashed",
                         lw=lw, c="k")
            else
                plt.text(x_right + 0.62, t, "\$t + E\$", fontsize=10,
                         ha="center", va="center")
                plt.plot([x_right + 0.94, 3], [t, t], linestyle="dashed",
                         lw=lw, c="k")
            end
        end
    end

    # time label
    plt.text(2.84, 2.47, "\$t\$", fontsize=10)

    # format
    ax.invert_yaxis()
    plt.yticks([0, 1, 2])
    plt.xticks([])
    plt.ylim([2.17, -0.22])
    plt.xlim([minimum(values(x_locs)) - 0.32, maximum(values(x_locs)) + 0.35])
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    for side in ["bottom", "top", "left"]
        ax.spines[side].set_color("#FFFFFF00")
    end
    ax.spines["right"].set_bounds([2.4, 0])
    ax.set_aspect("equal")
    plt.subplots_adjust(left=0.07, bottom=0.1, right=0.9, top=1.0)
    return (fig, ax)
end

function update_partitions(partitions, tree)
    p = partitions[end]
    info = get_tree_info(p["tree"])
    cells = MondrianForests.get_cells(tree)
    leaves = MondrianForests.get_cells(p["tree"])
    current_split = !isnothing(p["current"]) && !(p["current"] in cells)
    current_parent = !(p["current"] in leaves)

    # update time and tree
    if current_split && !current_parent
        new_time = minimum(t for t in times if t > p["time"])
        new_tree = MondrianForests.restrict(tree, new_time)
    else
        new_time = p["time"]
        new_tree = p["tree"]
    end

    # update current
    if isnothing(p["current"]) || current_parent || !current_split
        ids = [i[4]
               for i in info
               if i[3] in leaves &&
                  !(i[3] in p["terminals"]) && !(i[3] == p["current"])]
        ids = [i for i in ids if length(i) == minimum(length(j) for j in ids)]
        if !isempty(ids)
            new_current_id = minimum(ids)
            new_current = [i[3] for i in info if i[4] == new_current_id][]
        else
            new_current = nothing
        end
    else
        new_current = p["current"]
    end

    # update terminals
    if current_split || isnothing(p["current"])
        new_terminals = p["terminals"]
    else
        new_terminals = [p["terminals"]; [p["current"]]]
    end

    new_partition = Dict("time" => new_time,
                         "tree" => new_tree,
                         "current" => new_current,
                         "terminals" => new_terminals)

    return [partitions; [new_partition]]
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

# get locations of tree nodes for diagram
info = get_tree_info(tree)
ids = [i[1] for i in info]
xs = get_horizontal_value.(ids)
xs = invperm(sortperm(xs)) / 3
x_locs = Dict()
for i in 1:length(ids)
    x_locs[ids[i]] = xs[i]
end

# calculate the current and terminal nodes
partitions = [Dict("time" => 0.0,
                   "tree" => MondrianForests.restrict(tree, 0.0),
                   "current" => nothing,
                   "terminals" => [])]

info = get_tree_info(tree)
times = [i[2] for i in info]
for rep in 1:11
    global partitions = update_partitions(partitions, tree)
end

# plot the tree structures
println("plotting trees")
dpi = 500
for i in 1:length(partitions)
    println(i)
    partition = partitions[i]
    global (fig, ax) = plot_mondrian_tree(partition)
    plt.savefig("replication/construction_diagrams/mondrian_tree_$(i).png", dpi=dpi)
    plt.close("all")
end

# plot the generation of the partition
println("plotting partitions")
for i in 1:length(partitions)
    println(i)
    partition = partitions[i]
    global (fig, ax) = plot_mondrian_process(partition)
    plt.savefig("replication/construction_diagrams/mondrian_process_$(i).png", dpi=dpi)
    plt.close("all")
end
