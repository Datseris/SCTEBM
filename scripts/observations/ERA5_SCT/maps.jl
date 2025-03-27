# Visualize data as spatial maps
include("setup.jl")
# how to aggregate time dimension
# f(X) = timemean(X)
# date = "$(year(first(dims(C, Time)))) - $(year(last(dims(C, Time))))"

ti = 9
date = dims(C, Time)[ti]
f(X) = gnv(X[Time(ti)])
c = f(C)

figin = Figure(size = (1200, 1200))
figout = Figure(size = (1200, 1200))
# 3 by 3 layout

for k in 1:K # max number of in or out fields fields
    i, j = Tuple(CartesianIndices((3,5))[k])
    for (data, fig) in zip((inputs, outputs), (figin, figout))
        k > length(data) && continue
        X = data[k]
        x = f(X)
        x[oidxs] .= NaN
        ax, hm = heatmap(fig[i, j], x; colormap = Reverse(:dense))
        n = keys(data)[k]
        ax.title = String(n)
        ax.titlefont = :regular
        hideydecorations!(ax)
        ylims!(ax, 0, nothing)
    end
end

figuretitle!(figin, "ERA5 - inputs - $(nameprefix) - $(date)")
figuretitle!(figout, "ERA5 - outputs - $(nameprefix) - $(date)")

display(figin)
sleep(0.1)
display(figout)