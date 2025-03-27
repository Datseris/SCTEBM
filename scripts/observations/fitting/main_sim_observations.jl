# Script for running a CTMLM simulation given some input parameter values
# as well as specific parameterizations.
# It outputs a bunch of information regarding the system attractors,
# basins, global stability, as well as recording a bunch of observables.

# TODO: Make this script educative by expanding here the `process_multiparam_multistability_analysis` function!!!

using DrWatson
@quickactivate
using DynamicalSystems
using ConceptualClimateModels
using Statistics
using Distributions
include(srcdir("ctmlm_setups.jl"))
include(srcdir("simulations_run.jl"))
include(srcdir("simulations_process.jl"))

###########################################################################################
# Input
###########################################################################################
# Stage 1: load observational data as input parameters
# this ensures we use exactly the same combinations
include(scriptsdir("observations", "ERA5_SCT", "generate_observations.jl"))
# and this ensures that we discard points that are invalid
inputs = Dict(k => v[invalid_idxs] for (k, v) in inputs)

# Translate some of the observation quantities into the model parameters
D = inputs[:D] .* 1e-6
U = inputs[:U]
RH₊ = inputs[:RH₊]
δ_Δ₊T = inputs[:EIS]
S_0 = inputs[:S] # insolation
T_FTR = inputs[:T_FTR]
δ_FTR = (maximum(T_FTR) .- T_FTR)/2

input_fields = @ntuple D U RH₊ δ_FTR δ_Δ₊T S_0 α_a c_d

# and finally transform it to the specific form of the parameter container
# (vector of dictionaries)
param_container = map(i -> Dict(n => v[i] for (n, v) in pairs(input_fields)), eachindex(input_fields[1]))

# Stage 2: equations configuration
# This also defines the "default" model configuration
cooling = :q_x
cdversion = :sigmoid
invfix = :difference
Ld = :three_layer # has massive impact on the dynamics
ΔF = :three_layer # has massive impact on the dynamics
invdec = true
entrain = :Stevens2006
ftrgrad = :weak
input = @dict(cooling, cdversion, invfix, invdec, Ld, ΔF, entrain, ftrgrad)

# Stage 3: folder name to save to
foldergroup = ["sims", "observed_params"]

# The rest is the automated running code
# Keep in mind that parameter fitting is done inside the `ctmlm_setup` function.
process_multiparam_multistability_analysis(
    input, param_container, foldergroup; force = true,
    kw_analysis = (density = 5, extra = 2, N = 1000),
    kw_process = (plot_rmi = false, delete_files = true, estimate_rmi = false),
)
