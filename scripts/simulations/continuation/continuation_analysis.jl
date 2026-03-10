# This script loads a continuation (climate change scenario simulations)
# and analyzes it by plotting any user-chosen set of variables or observables
# of the dynamical system. Simply change the `observables` variable below
# to plot other observables. If the simulation does not exist (combination of `input`
# and `used`, just run the `continuation_generation.jl` script.

using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
using Statistics
using Distributions
using Statistics: mean
include(srcdir("sctebm_setups.jl"))
include(srcdir("simulations_process.jl"))
include(srcdir("theme.jl"))
include("continuation_plotting_helpers.jl")
TRANSPARENCY[] = 0.25

# Inputs to visualize:
cooling = :q_x
invfix = :temperature
Ld = :three_layer
ΔF_s = :ctrc
co2 = 2
cloud_rad = :lwp
entrain = :Stevens2006
ftrgrad = :none
input = @dict(cooling, invfix, Ld, ΔF_s, co2, entrain, ftrgrad, cloud_rad)

used = sort([:D, :CO2]) # always sort this!!!

foldergroup = ["sims", "continuations"]
prefix = "used="*join(string.(used), "+")
name = savename(input)*"_"*prefix
data = wload(datadir(foldergroup..., name)*".jld2")
@unpack continuations, param_values = data

# observables = [:LHF, :RCT, :CTRC, :Λ, :C, :z_b]
observables = [:SST, :C]

ids = 1:10

# plotsy plotsy
fig = Figure(size = (figwidth/2, 2figheight))
axs = load_n_plot_continuation!(fig, input, used, observables; ids)
# ylims!(axs[4], 0, 8)
resize!(fig, 600, 400)
display(fig)

# wsave(papersdir("figures", "continuation_analysis"), fig)
# %%
fig = Figure(size = (figwidth/2, figheight))
ax = Axis(fig[1,1])
load_n_plot_continuation_only_c!(ax, input, used)
figuretitle!(fig, string(used))
fig