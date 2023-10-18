using PyPlot

# TODO makefile to minify then copy to docs

# plot setup
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["text.latex.preamble"] = "\\usepackage{montserrat}"
plt.ioff()

# start plot
(fig, ax) = plt.subplots(figsize=(3, 3))
for side in ["bottom", "top", "left", "right"]
    ax.spines[side].set_color("#FFFFFF00")
end
plt.xticks([])
plt.yticks([])

# curve
n = 200
xs = range(0, 1.005, length=n)
ys = (1.2 .* xs .- 0.63).^3 .+ 0.55 - 0.2 * xs.^2 -
0.05 * (0.5 .* xs .+ 0.5).^10 +
0.07 * (1 .- xs).^10
plt.plot(xs, ys, c="k", lw=2)

# box
plt.plot([0, 1], [0, 0], c="k", lw=2)
plt.plot([0, 0], [0, ys[1] - 0.005], c="k", lw=2)
plt.plot([1, 1], [0, ys[end] - 0.006], c="k", lw=2)

# text
plt.text(0.49, 0.59, "\\bfseries Mondrian Forests", fontsize=19, ha="center", va="bottom")

# piet lines
plt.plot([0.25, 0.25], [0, 0.498], c="k", lw=2)
plt.plot([0, 1], [0.25, 0.25], c="k", lw=2)
plt.plot([0.85, 0.85], [0, 0.25], c="k", lw=2)
plt.plot([0.85, 1.0], [0.1, 0.1], c="k", lw=2)

# piet blocks
red_xs = [xs[i] for i in 1:n if xs[i] >= 0.25]
red_ys = [ys[i] for i in 1:n if xs[i] >= 0.25]
plt.fill_between(red_xs, 0.25, red_ys, fc="#FA0000")
plt.fill_between([0, 0.25], 0, 0.25, fc="#0E75CF")
plt.fill_between([0.85, 1], 0, 0.1, fc="#FFEC02")

# trim
plt.plot([-0.011, -0.011], [0, 0.55], c="w", lw=2)
plt.plot([1.011, 1.011], [0, 0.65], c="w", lw=2)

# save
fig.subplots_adjust(bottom=0.08, top=0.9, left=0.07, right=0.93)
plt.savefig("./replication/logo/logo.svg", dpi=1000)
plt.close("all")
