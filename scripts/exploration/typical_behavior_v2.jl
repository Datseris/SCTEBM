using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
using Statistics
using Distributions
using Statistics: mean
include(srcdir("sctebm_setups.jl"))
include(srcdir("simulations_run.jl"))
include(srcdir("theme.jl"))

# create random initial conditions
Cgrid = [0.05, 0.5, 0.95]
Tgrid = [285.0, 290.0, 295.0, 300.0]
u0s = random_initial_conditions(Cgrid, Tgrid, 5)

colors = color_from_u0.(u0s)

function overplot!(ds, ax, limitcycle = false) # when to plot endpoints
    for (i, u0) in enumerate(u0s)
        X, t = trajectory(ds, 150.0, u0; Δt = 0.05, save_idxs = [:C, :SST])
        # it looks worse with the initial conditions
        # scatter!(ax, X[1]; color = :white, strokewidth = 1, markersize = 5, strokecolor = :black)
        fadelines!(ax, X; color = colors[i], linewidth = 2, fade = 0.3, linestyle = :solid)
        if !limitcycle # plot all endpoints
            scatter!(ax, X[end]; marker = :circle, color = :red, strokewidth = 2, markersize = 15, strokecolor = :black)
        end
        if i == length(u0s) && limitcycle
            X, t = trajectory(ds, 25.0, u0; Δt = 0.55, Ttr = 300, save_idxs = [:C, :SST])
            scatterlines!(ax, X; color = :red, strokewidth = 2, markersize = 15, strokecolor = :black)
        end
    end
end

# %% setup simulation
# Start at a bistable configuration
ds, eqs = sctebm_setup(;
    cooling = :q_x,
    cdversion = :sigmoid,
    invfix = :difference,
    Ld = :three_layer,
    ΔF_s = :three_layer,
    invdec = true,
    entrain = :Stevens2006,
    ftrgrad = :weak,
    cloud_rad = :lwp,
    cloud_limit = true,
)

set_parameter!(ds, :U, 5.3)
set_parameter!(ds, :D, 3e-6)
set_parameter!(ds, :α_C_max, 0.8)
set_parameter!(ds, :δ_Δ₊T, 7.0)

# hard to find parameters that do _everything_ for exactly the same cloud setup...
params = [
    [:D => 3.5e-6, :U => 3.5, :δ_FTR => 10, :δ_Δ₊T => 10, :RH₊ => 0.35, :α_C_max => 0.7],
    [:D => 3e-6, :U => 5.3, :δ_FTR => 0.0, :δ_Δ₊T => 7, :RH₊ => 0.2, :α_C_max => 0.8],
    [:D => 2e-6, :U => 9.0,], # Cu
    [:D => 5.5e-6, :U => 4.0,], # cloud-free
    [:D => 3.88e-6, :U => 5.56, :δ_Δ₊T => 4.8, :RH₊ => 0.1],
]

titles = ["Sc only", "Bistable", "Cu only", "Cloud-free", "Tristable"]

# plot
fig, axs = axesgrid(2, 3; xlabels = "C", ylabels = "SST [K]", sharey = true, sharex = true)

for (j, p) in enumerate(params)
    ax = transpose(axs)[j]
    set_parameters!(ds, p)
    overplot!(ds, ax)
    ax.title = titles[j]
    ylims!(ax, 285, 305)
end

# add the limit cycle
ds, eqs = sctebm_setup(;
    Ld = :three_layer,
    ΔF_s = :ctrc,
    entrain = :Stevens2006,
    ftrgrad = :weak,
    cloud_rad = :const,
    cloud_limit = true,
)

set_parameter!(ds, :U, 9)
set_parameter!(ds, :D, 4.5e-6)
set_parameter!(ds, :RH₊, 0.2)
set_parameter!(ds, :δ_Δ₊T, 5)
set_parameter!(ds, :δ_FTR, 5)

overplot!(ds, axs[end], true)
axs[end].title = "Limit cycle (unphysical)"
ylims!(axs[end], 285, 305)

display(fig)

# %% final decorations
label_axes!(axs; transpose = true)
xs = -0.05:0.005:0.05
ys = -7.5 .* xs.^2
streamlines!(fig.scene, xs .+ 0.36, ys .+ 0.96; space = :relative, linestyle = :dash, color = "black", markersize = 20)
text!(fig.scene, 0.36 + maximum(xs), 0.95; text = "± envir.", space = :relative, )
streamlines!(fig.scene, xs .+ 0.68, ys .+ 0.96; space = :relative, linestyle = :dash, color = "black", markersize = 20)
text!(fig.scene, maximum(xs) + 0.68, 0.95; text = "± envir.", space = :relative, )

display(fig)

# wsave(papersdir("figures", "typical_behavior.png"), fig)
