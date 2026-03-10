# Quick and simple script that given a SCTEBM setup it produces distribution of
# some chosen observables for some distribution of input parameters.
using DrWatson
@quickactivate
using DynamicalSystems
using Distributions
include(srcdir("sctebm_setups.jl"))
include(srcdir("simulations_process.jl"))

input = Dict(
    :cooling => :q_x,
    :entrain => :Stevens2006,
    :cloud_rad => :lwp,
    :cloud_limit => true,
    :invfix => :difference,
    :ftrgrad => :weak,
    :Ld => :three_layer,
    :ΔF_s => :three_layer,
)

distributions = Dict(
    :D => Uniform(1e-6, 5e-6),
    :U => Uniform(5.0, 9.0), # the exaggerated range of U takes into account d variability
    :δ_Δ₊T => Uniform(1, 10),
    :RH₊ => Uniform(0.1, 0.4),
    :δ_FTR => Uniform(0, 10),
)

observables_to_obtain = [:C, :SST, :α_C, :ε_C, :LWP]

# Plot distributions
foldergroup = ["sims", "quick_obs_distr"]
output = process_multiparam_multistability_analysis(
    input, distributions, foldergroup;
    observables_to_obtain, force = true,
    # multistability multiparameter keywords
    kw_analysis = (density = 3, extra = 2, N = 200),
    # keywords for processing
    kw_process = (estimate_rmi = false, plot_density = true, save_plots = false, delete_files = true, save_processed_output = false),
)
