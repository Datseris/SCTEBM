# Same as `playground_setups.jl` but starting already at a limit cycle.

using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
include(srcdir("ctmlm_setups.jl"))
include("playground_helpers.jl")

ds, eqs = ctmlm_setup(;
    cooling = :q_x,
    Ld = :three_layer,
    ΔF = :three_layer,
    invfix = :difference,
    ftrgrad = :none,
    entrain = :Stevens2006,
    starting_parameters = Dict(:D => 2.5e-6, :U => 7, :q_x_rate => 2, :RH₊ => 0.1)
)

# %% Configure and launch GUI
GUI_obs = [:SST, :q_b, :s_b, :z_b, :T₊]
GUI_par = [:U, :D, :δ_Δ₊T, :δ_FTR, :q_x_rate, :RH₊]
fig, dsobs = ctmlm_gui(ds, GUI_par, GUI_obs)

# %%
figuretitle!(fig, "SCTEBM exploration via a GUI")
wsave(plotsdir("limitcycle_gui"), fig)