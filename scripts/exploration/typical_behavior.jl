using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
using Statistics
using Distributions
using Statistics: mean
include(srcdir("CloudToppedBulkBoundaryLayer", "CloudToppedBulkBoundaryLayer.jl"))
include(srcdir("ctbbl_setups.jl"))
include(srcdir("theme.jl"))
using MakieExtra

using Random
Random.seed!(1234)
# create random initial conditions
u0s0 = [Dict(:C => C, :SST => T) for C in [0.05, 0.5, 0.95] for T in [285.0, 290.0, 295.0, 300.0]]
u0s = Dict{Symbol, Real}[]
for u0 in u0s0
    for _ in 1:5
        # get random variability coefficientt
        r = () -> 1 + 1e-2randn()
        T0 = u0[:SST]
        RH0 = clamp(0.7r(), 0, 1)
        u0[:SST] = T0*(0.5 + 0.5r())
        u0[:C] = u0[:C]*r()
        u0[:s_b] = (T0 - 2.0)*r() # in equilibrium s is about 1-2 Kelvin less than SST.
        u0[:q_b] = RH0*CTMLM.q_saturation(T0)
        u0[:z_b] = 1200*r()
        push!(u0s, copy(u0))
    end
end

colors = color_from_u0.(u0s)

function overplot!(ds, ax, when = true) # when to plot endpoints
    for (i, u0) in enumerate(u0s)
        X, t = trajectory(ds, 100.0, u0; Δt = 0.05, save_idxs = [:C, :SST])
        # it looks worse with the initial conditions
        # scatter!(ax, X[1]; color = :white, strokewidth = 1, markersize = 5, strokecolor = :black)
        fadelines!(ax, X[:, 1], X[:, 2]; color = colors[i], linewidth = 2, fade = 0.3)
        if when # plot all endpoints
            scatter!(ax, X[end]; color = :red, strokewidth = 2, markersize = 15, strokecolor = :black)
        end
        if i == length(u0s) && !when
            X, t = trajectory(ds, 25.0, u0; Δt = 0.55, Ttr = 300, save_idxs = [:C, :SST])
            scatterlines!(ax, X; color = :red, strokewidth = 2, markersize = 15, strokecolor = :black)
        end
    end
end

# setup simulation
cooling = :q_x
invfix = :difference
Ld = :three_layer
entrain = :Stevens2006
ΔF_s = :crc
input = @dict(cooling, invfix, w_m, Ld, ΔF_s, entrain)

ds, eqs = ctbbl_setup(; input...)
set_parameter!(ds, :RH₊, 0.1)
set_parameter!(ds, :D, 3e-6)
params = [5.5, 8.0, 10.5]

# plot
fig, axs = axesgrid(2, 2; size = (600, 600), xlabels = "C", ylabels = "SST", sharey = true, sharex = true)

for (j, p, ax) in zip(1:3, params, axs)
    set_parameter!(ds, :U, p)
    reinit!(ds)
    overplot!(ds, ax, j != 2)
    ylims!(ax, 285, 305)
end

# %% plot only bistable regime
empty!(axs[4])
input[:ΔF_s] = :three_layer
input[:entrain] = :Stevens2006
input[:ftrgrad] = :weak
input[:Ld] = :three_layer
ds, eqs = ctbbl_setup(; input...)
set_parameter!(ds, :U, 4.5)
set_parameter!(ds, :D, 4e-6)
set_parameter!(ds, :RH₊, 0.2)
set_parameter!(ds, :SST_X_0, 0)
set_parameter!(ds, :α_c, 0.4)
set_parameter!(ds, :q_x_rate, 2.0) # lowering this was the trick!
reinit!(ds)
overplot!(ds, axs[4], true)
ylims!(axs[4], 285, 305)
fig |> display

# %%
# arrows over whole figure
ax = Axis(fig[:, :])
ax.limits = ((-2, 2), (0, 2))
xmin = -0.6
xs = xmin:0.01:(-xmin)
ys = xs.^2 .+ 0.8

arrowlines!(ax, xs, ys; arrowstyle = "|>-|>", markersize = 15, color = :black)
text!(ax, xs[end] .+ 0.05, ys[end]; text = "changing\nenvironmental\nconditions", align = (:left, :center))

hidedecorations!(ax)

hidespines!(ax)

label_axes!(axs)
fig |> display

wsave(papersdir("figures", "typical_behavior.png"), fig)
