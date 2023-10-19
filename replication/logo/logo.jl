using PyPlot

# start plot
plt.ioff()
(fig, ax) = plt.subplots(figsize=(5, 5))
for side in ["bottom", "top", "left", "right"]
    ax.spines[side].set_color("#FFFFFF00")
end
plt.xticks([])
plt.yticks([])

# colors
col_dark = "#111111"
col_light = "#f9f9ee"
col_blue = "#4063d8"
#col_green = "#389826"
col_green = "#38a826"
#col_purple = "#9558b2"
col_purple = "#a588e2"
col_red = "#cb3c33"

# outer box
lw = 12
x_eps = 0.037
y_eps = 0.020
top = 0.59
xs_outer = [-x_eps, 1+x_eps, 1+x_eps, -x_eps, -x_eps]
ys_outer = [-y_eps, -y_eps, top+y_eps, top+y_eps, -y_eps]
outer_col = col_purple
plt.plot(xs_outer, ys_outer, c=outer_col, lw=2*lw, zorder=0)
plt.fill(xs_outer, ys_outer, c=outer_col, lw=0)

# curve
n = 100
xs = range(0, 1, length=n)
ys = (1.2 .* xs .- 0.63).^3 .+ 0.60 - 0.2 * xs.^2 -
0.05 * (0.5 .* xs .+ 0.5).^10 +
0.07 * (1 .- xs).^10
all_xs = [[0.5, 0]; xs; [1, 0.5]]
all_ys = [[0, 0]; ys; [0, 0]]
plt.plot(all_xs, all_ys, c="#111111", lw=lw)

# piet lines
eps = 0.01
x1 = 0.25
x2 = 0.80
y1 = 0.11
y2 = 0.25
plt.plot([x1, x1], [eps, 0.55], c="#111111", lw=lw)
plt.plot([eps, 1-eps], [y2, y2], c="#111111", lw=lw)
plt.plot([x2, x2], [eps, y2-eps], c="#111111", lw=lw)
plt.plot([x2+eps, 1-eps], [y1, y1], c="#111111", lw=lw)

# piet blocks
red_xs = [xs[i] for i in 1:n if xs[i] >= x1]
red_ys = [ys[i] for i in 1:n if xs[i] >= x1]
plt.fill_between(red_xs, y2, red_ys, fc=col_blue)
white_xs = [xs[i] for i in 1:n if xs[i] <= x1]
white_ys = [ys[i] for i in 1:n if xs[i] <= x1]
plt.fill_between(white_xs, y2, white_ys, fc=col_light)
plt.fill_between([0, x1], 0, y2, fc=col_red)
plt.fill_between([x2, 1], 0, y1, fc=col_green)
plt.fill_between([x1, x2], 0, y2, fc=col_light)
plt.fill_between([x2, 1], y1, y2, fc=col_light)

# save
plt.savefig("./docs/src/assets/logo.svg", dpi=1000, transparent=true)
plt.close("all")
