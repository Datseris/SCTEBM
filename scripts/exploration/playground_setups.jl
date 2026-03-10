# Script for playing around with different equations and parameters for CTMLM
# There are two ways to specify the model equations;
# 1. use the `sctebm_setup` function (this file) that makes same setups as in the paper
# 2. directly create your own equations (utilizing the default equations in defaults.jl),
# which is in the file `playground_equations.jl`

# From there, you specify what observables to visualize, and for which
# parameters to create sliders for. Both of these options are given
# after model initialization, which allows one to change change parameter
# sliders or observables visualized without restarting the model!

using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
include(srcdir("sctebm_setups.jl"))
include("playground_helpers.jl")

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

# %% Configure and launch GUI
# the `GUI_obs` can be a premade configuration or a vector of things to observe
# which can include arbitrary symbolic expressions or functions as per DynamicalSystems.jl.
# See `playground_helpers.jl` for predefined options
GUI_obs = :none
# GUI_obs = [:CTRC, :LHF, :RCT, :Λ, :α_C, :ε_C]
# GUI_obs = [:CTRC, :SST, :α_C, :ε_C, :CLT, :z_b]
GUI_obs = [:CTRC, :α_C, :ε_C, :C_Λ, :C_κ, :CLT]

# The `GUI_par` can only be a vector of symbols each representing a system parameter
GUI_par = [:U, :D, :δ_Δ₊T, :g_C, :RH₊, :α_C_max]

fig, dsobs = ctmlm_gui(ds, GUI_par, GUI_obs)
