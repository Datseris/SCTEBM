# Script that generates the profiles versus height given an instance of the SCTEBM
using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
include(srcdir("sctebm_setups.jl"))
include(srcdir("simulations_run.jl"))
include(srcdir("theme.jl"))

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

@assert length(attractors) == 2

Cu = attractors[findfirst(A -> A[:, :C][end] < 0.5, attractors)]
Cu = Dict(s => Cu[1, s] for s in (:z_b, :s_b, :q_b, :C, :SST))
Sc = attractors[findfirst(A -> A[:, :C][end] > 0.5, attractors)]
Sc = Dict(s => Sc[1, s] for s in (:z_b, :s_b, :q_b, :C, :SST))

# Profiles
function profiles!(axs, ds, state::Dict)
    set_state!(ds, state)
    z_b, s_b, q_b, = observe_state.(ds, (:z_b, :s_b, :q_b, ))
    z_LCL, i_Λ, λ_q, λ_s = observe_state.(ds, (:z_cb, :i_Λ, :λ_q, :λ_s ))
    T₊, q₊, s₊, Δ₊q, Δ₊s = observe_state.(ds, (:T₊, :q₊, :s₊, :Δ₊q, :Δ₊s))

    q_LCL = i_Λ*λ_q*Δ₊q + q_b
    s_LCL = i_Λ*λ_s*Δ₊s + s_b

    T_0 = s_b
    T_LCL = CTMLM.temperature_exact(z_LCL, s_b, q_b)
    T_t = CTMLM.temperature_exact(z_b, s_b, q_b)

    ql_t = CTMLM.q_liquid(T_t, q_b, z_b)
    # profiles
    zgrid = [0, z_LCL, z_LCL, z_b, z_b, 1.25z_b]
    lines!(axs[1], [s_b, s_b, s_LCL, s_LCL, s₊, s₊], zgrid)
    lines!(axs[2], [q_b, q_b, q_LCL, q_LCL, q₊, q₊], zgrid; color = Cycled(2))
    lines!(axs[3], [0, 0, 0, ql_t, 0, 0], zgrid; color = Cycled(3))
    lines!(axs[4], [T_0, T_LCL, T_LCL, T_t, T₊, T₊], zgrid; color = Cycled(4))
    ylims!(axs[4], 0, nothing)
    axs[1].yticklabelrotation = π/2
    for ax in axs
        ax.yticks = ([0, z_LCL, z_b], ["0", "$(round(Int, z_LCL))", "$(round(Int, z_b))"])
    end
    return
end

fig, axs = axesgrid(2, 4; xlabels = ["s [K]", "q [g/kg]", "q_l [g/kg]", "T [K]"], ylabels = ["z [m], Sc state", "z [m], Cu state"], sharey = true, sharex = true, size = (figwidth, 2figheight))

profiles!(axs[1, 1:4], ds, Sc)
profiles!(axs[2, 1:4], ds, Cu)
display(fig)
wsave(papersdir("figures", "profiles"), fig)
