# This script explores the connection between liquid water potential
# temperature and the decoupling variable Λ!

using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
include(srcdir("sctebm_setups.jl"))
include(srcdir("simulations_run.jl"))

# Make a bistable configuration
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
set_parameter!(ds, :D, 3e-6) # make this 4 for tristable
set_parameter!(ds, :α_C_max, 0.8)
set_parameter!(ds, :δ_Δ₊T, 7.0)

# Get attractors (steady states)
fractions, labels, convergence, attractors = multistability_analysis(ds)

Cu = attractors[findfirst(A -> A[:, :C][end] < 0.5, attractors)]
Cu = Dict(s => Cu[1, s] for s in (:z_b, :s_b, :q_b, :C, :SST))
Sc = attractors[findfirst(A -> A[:, :C][end] > 0.5, attractors)]
Sc = Dict(s => Sc[1, s] for s in (:z_b, :s_b, :q_b, :C, :SST))

function theta_liquid(s, q, z)
    T = CTMLM.temperature_exact(z, s, q)
    p_z = CTMLM.pressure(z, T)
    Π = (p_z/CTMLM.p₀)^0.286
    θl = (s - CTMLM.g*z/CTMLM.cₚ)/Π
    return θl
end

function theta_ratio(s_b, q_b)
    θl_0 = theta_liquid(s_b, q_b, 0)
    z_LCL = CTMLM.cloud_base_height_bolton1980(s_b, q_b)
    θl_LCL = theta_liquid(s_b, q_b, z_LCL)
    return (θl_LCL - θl_0)/θl_0
end

function theta_ratio_alt(ds, state::Dict)
    @unpack q_b, s_b = state
    θl_0 = theta_liquid(s_b, q_b, 0)
    set_state!(ds, state)
    z_LCL, i_Λ, λ_q, λ_s, Δ₊q, Δ₊s = observe_state.(ds, (:z_cb, :i_Λ, :λ_q, :λ_s, :Δ₊q, :Δ₊s ))
    q_LCL = i_Λ*λ_q*Δ₊q + q_b
    s_LCL = i_Λ*λ_s*Δ₊s + s_b
    θl_LCL = theta_liquid(s_LCL, q_LCL, z_LCL)
    return (θl_LCL - θl_0)/θl_0
end

# %% plot
fig, axs = axesgrid(2,2; titles = ["Sc", "Cu"], xlabels = "decoupling Λ", ylabels = ["θl ratio", "θl ratio - alt."], sharex = true, size = (600, 600))

for (i, state) in enumerate((Sc, Cu))
    set_state!(ds, state)
    Λ = observe_state(ds, :Λ)
    θl = theta_ratio(state[:s_b], state[:q_b])
    θl2 = theta_ratio_alt(ds, state)
    scatter!(axs[1, i], Λ, θl)
    scatter!(axs[2, i], Λ, θl2)
end

fig

