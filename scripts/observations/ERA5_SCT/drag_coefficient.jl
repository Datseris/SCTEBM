include("setup.jl")
using LsqFit

V = inputs[:U]
d = inputs[:drag]*1000

# %%
fig, axs = axesgrid(1, 2)

bad1 = findall(x -> ismissing(x) || isnan(x) || isinf(x), vec(d))
bad2 = findall(x -> ismissing(x) || isnan(x) || isinf(x), vec(V))
bad = sort!(union(bad1, bad2))

y = deleteat!(copy(vec(d)), bad)
x = deleteat!(copy(vec(V)), bad)

density!(axs[1], y)
axs[1].ylabel = "density"
axs[1].xlabel = "drag*1000"

hideydecorations!(axs[1]; label = false)
scatter!(axs[2], x, y)
axs[2].ylabel = "drag*1000"
axs[2].xlabel = "V"

display(fig)
