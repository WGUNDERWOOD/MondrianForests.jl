using CSV
using DataFrames
using PyPlot
using Random
using Colors
using Plots
using Dates

Random.seed!(314159)
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
# plt.rc('text', usetex=True)
rcParams["text.latex.preamble"]="\\usepackage[sfdefault,light]{FiraSans}"
plt.ioff()

function normalize_data(x)
    eps = 0.001
    x_new = (x .- minimum(x)) / (maximum(x) - minimum(x))
    x_new = (1 - 2 * eps) .* x_new .+ eps
end

# load data
data = DataFrame(CSV.File("replication/weather.csv", limit=nothing, missingstring="NA"))
#data = data[Date("2016-01-01") .<= data.Date .< Date("2017-01-01"), :]
data = data[:, [:Pressure3pm, :Humidity3pm, :RainTomorrow]]
dropmissing!(data)
data.RainTomorrow = replace(data.RainTomorrow, "Yes" => 1, "No" => 0)
data = Float64.(data)
println(sum(data.RainTomorrow .== 0))
println(sum(data.RainTomorrow .== 1))

# process data
pres_lims = [995, 1035]
data = data[pres_lims[1] .<= data.Pressure3pm .<= pres_lims[2], :]
data.Humidity3pm .+= 1 * (rand(nrow(data)) .- 0.5)
data.Pressure3pm .+= 0.5 * (rand(nrow(data)) .- 0.5)

# plot data
(fig, ax) = plt.subplots(figsize=(4, 3.5))
dry_color = "#f54e00"
wet_color = "#0045dd"
plt.scatter(data.Humidity3pm, data.Pressure3pm,
            c = [[dry_color, wet_color][Int(data.RainTomorrow[i]) + 1]
                 for i in 1:nrow(data)],
            s=5, alpha=0.4, marker=".", ec=nothing)

# format plot
#plt.xticks([0, 1], labels=hum_lims)
plt.yticks([1000 + i * 10 for i in 0:3])
plt.xlabel("Relative humidity at 3pm (\\%)")
plt.ylabel("Pressure at 3pm (mbar)")
#ax.xaxis.set_label_coords(0.5, -0.04)
#ax.yaxis.set_label_coords(-0.04, 0.5)
#ax.tick_params(color="w", direction="in", pad=0)
#for side in ["bottom", "top", "left", "right"]
    #ax.spines[side].set_color("#FFFFFF00")
#end

# color key
dry_handle = plt.scatter([], [], c=dry_color)
wet_handle = plt.scatter([], [], c=wet_color)
ax.legend([dry_handle, wet_handle],
          ["Dry tomorrow", "Wet tomorrow"],
          frameon=false, bbox_to_anchor=(0.95, 1.16), ncol=2)

plt.tight_layout()
#PyPlot.savefig("replication/weather_data.pgf")
PyPlot.savefig("replication/weather_data.pdf")
PyPlot.savefig("replication/weather_data.png", dpi=300)
plt.close("all")
