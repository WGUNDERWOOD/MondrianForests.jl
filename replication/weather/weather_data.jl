using CSV
using DataFrames
using PyPlot
using Random
using Colors
using Plots
using Dates
using MondrianForests

# plot setup
Random.seed!(314159)
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["text.latex.preamble"]="\\usepackage[sfdefault,light]{FiraSans}"
plt.ioff()

function load_data(limit=nothing)
    file = "replication/weather/weather.csv"
    data = DataFrame(CSV.File(file, limit=limit, missingstring="NA"))
    data = data[:, [:Humidity3pm, :Pressure3pm, :RainTomorrow]]
    dropmissing!(data)
    data.RainTomorrow = replace(data.RainTomorrow, "Yes" => 1, "No" => 0)
    data = Float64.(data)
    pres_lims = [985, 1040]
    data = data[pres_lims[1] .<= data.Pressure3pm .<= pres_lims[2], :]
    data.Humidity3pm .+= 1 * (rand(nrow(data)) .- 0.5)
    data.Pressure3pm .+= 0.5 * (rand(nrow(data)) .- 0.5)
    return data
end

function normalize(x)
    return (x .- minimum(x)) / (maximum(x) - minimum(x))
end

function plot_data(data)

    # main plot
    (fig, ax) = plt.subplots(figsize=(4, 3.5))
    dry_color = "#d55e00"
    wet_color = "#0035dd"
    colors = [dry_color, wet_color][Int.(data.RainTomorrow) .+ 1]
    plt.scatter(data.Humidity3pm, data.Pressure3pm, c=colors,
                s=5, alpha=0.3, marker=".", ec=nothing)
    plt.yticks([990 + i * 10 for i in 0:5])
    plt.xlabel("Relative humidity at 3pm (\\%)")
    plt.ylabel("Pressure at 3pm (mbar)")

    # color key
    dry_handle = plt.scatter([], [], c=dry_color)
    wet_handle = plt.scatter([], [], c=wet_color)
    ax.legend([dry_handle, wet_handle], ["Dry tomorrow", "Wet tomorrow"],
              handletextpad=0.1, frameon=false,
              bbox_to_anchor=(0.93, 1.16), ncol=2)

    plt.tight_layout()
    PyPlot.savefig("replication/weather_data.png", dpi=300)
    plt.close("all")
end

function plot_data_partition(data, tree)
    splits = MondrianForests.get_splits(tree)
    cells = MondrianForests.get_cells(tree)
    (fig, ax) = plt.subplots(figsize=(4, 3.5))

    # plot root cell
    lw = 0.3
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

    # main plot
    (fig, ax) = plt.subplots(figsize=(4, 3.5))
    dry_color = "#d55e00"
    wet_color = "#0035dd"
    colors = [dry_color, wet_color][Int.(data.RainTomorrow) .+ 1]
    plt.scatter(data.Humidity3pm, data.Pressure3pm, c=colors,
                s=5, alpha=0.3, marker=".", ec=nothing)
    plt.yticks([990 + i * 10 for i in 0:5])
    plt.xlabel("Relative humidity at 3pm (\\%)")
    plt.ylabel("Pressure at 3pm (mbar)")

    # color key
    dry_handle = plt.scatter([], [], c=dry_color)
    wet_handle = plt.scatter([], [], c=wet_color)
    ax.legend([dry_handle, wet_handle], ["Dry tomorrow", "Wet tomorrow"],
              handletextpad=0.1, frameon=false,
              bbox_to_anchor=(0.93, 1.16), ncol=2)

    plt.tight_layout()
    PyPlot.savefig("replication/weather_data_partition.png", dpi=300)
    plt.close("all")
end

data = load_data()
display(data)
#plot_data(data)
lambda = 3.0
tree = MondrianTree(2, lambda)
# TODO work out how to rescale into [0,1]^2
plot_data_partition(data, tree)

