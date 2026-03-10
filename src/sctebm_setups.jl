# Provide a central function `sctebm_setup` that creates
# a dynamical system representing the SCTEBM. Create any variant by utilizing
# keywords that map aspect to its options.

using ConceptualClimateModels
import ConceptualClimateModels.CloudToppedMixedLayerModel as CTMLM

"""
    sctebm_setup(; aspect = option, ...) → ds, eqs

Main function associated with the research article (Datseris 2026).
Given an input of SCTEBM aspect options, create and return the resulting dynamical system
and equations (that were given to `processes_to_coupledodes`).
The aspects are named slightly different from the Table 1 of the paper,
so it is best to inspect the source code of the function for their names.
"""
function sctebm_setup(;
        dt = 0.05, diffeq = (adaptive = false, dt, verbose = false),
        extra_eqs = [], starting_parameters::Dict = Dict(),
        # The keywords are the aspects and their values are the options:
        invfix = :difference, # :difference or :temperature
        ftrgrad = :weak, # :none, :weak:, :strong, or number
        Ld = :three_layer, # :three_layer, :fixed, or other options, see function
        ΔF_s = :ctrc, # :ctrc, :three_layer, :Gesso2014
        co2::Int = 1, # 1 2 or 3 for increasing CO2 impact
        entrain = :Gesso14, # :Gesso2014 or :Stevens2006
        cooling = :q_x,
        w_m = false,
        cdversion = :sigmoid,
        cloud_limit = true,
        cloud_rad = :const, # :const or :lwp
        invdec = true, # enable decoupling augmentations with iλ
        fcrc = true, # fractional cloud radiative effects, ∝ C
    )

    # Core equations
    eqs = [
        # Cloud equations
        CTMLM.cf_dynamic(cdversion; thinness_limiter = cloud_limit),
        CTMLM.decoupling_variable(),
        CTMLM.cloud_base_height(:Bolton1980),
        # Boundary layer
        CTMLM.sst_dynamic(),
        CTMLM.mlm_dynamic(),
        CTMLM.ρ₀ ~ 1.2, # just for simplicity and reducing compute
        CTMLM.mlm_q₊(:relative),
        CTMLM.mlm_s₊(invfix;
            cloud_effect = false, CO2_effect = co2 == 3,
        ),
        CTMLM.entrainment_velocity(entrain),
        ParameterProcess(CTMLM.s_x, 0.0),

        # Temperatures
        CTMLM.free_troposphere_emission_temperature(ftrgrad; add_co2 = true),
        CTMLM.cloud_emission_temperature(:top),
        CTMLM.bbl_emission_temperature(),

        # Radiation
        CTMLM.cloud_longwave_cooling(),
        CTMLM.cloud_shortwave_warming(:insolation; cloud_fraction = fcrc),
        CTMLM.mlm_radiative_cooling(ΔF_s),
        CTMLM.downwards_longwave_radiation(Ld),
    ]

    @parameters SST_X_0 = 0.0
    SST_X_rhs = SST_X_0
    if co2 > 1 # surface warming due to global warming
        # Half because it is unreasonable that all of TOA radiative forcing goes into SST
        SST_X_rhs += (3.6/2)*log2(CTMLM.CO2/400)
    end
    # Cooling:
    if cooling == :q_x
        @parameters q_x_rate = 1.5
        push!(eqs, CTMLM.q_x ~ q_x_rate*max(CTMLM.q_saturation(CTMLM.SST)/CTMLM.q_saturation(290) - 1, 0))
    elseif cooling == :sst_x
        @parameters SST_X_cooling = -2.0
        SST_X_rhs += SST_X_cooling*max(CTMLM.σ_SB*(CTMLM.SST^4 - 290^4), 0)
    end

    push!(eqs, CTMLM.SST_X ~ SST_X_rhs)

    if w_m
        @parameters e_m = 0.001
        push!(eqs, CTMLM.w_m ~ e_m*CTMLM.C)
    end

    if invdec
        append!(eqs, [
            # clamping makes numerical integration stable
            CTMLM.λ_q ~ clamp((CTMLM.CLT/2750)^1.3, 0, 0.5),
            CTMLM.λ_s ~ 0.5*CTMLM.λ_q,
            CTMLM.Δ₊q ~ (1 - CTMLM.i_Λ*CTMLM.λ_q)*(CTMLM.q₊ - CTMLM.q_b),
            CTMLM.Δ₊s ~ (1 - CTMLM.i_Λ*CTMLM.λ_s)*(CTMLM.s₊ - CTMLM.s_b),
        ])
    end

    # cloud radiation (always proportional to cloud fraction in the SCTEBM)
    if cloud_rad == :const
        append!(eqs, [
            CTMLM.cloud_emissivity(1.0),
            CTMLM.cloud_albedo(0.38),
        ])
    elseif cloud_rad == :lwp
        append!(eqs, [
            CTMLM.cloud_emissivity(:lwp),
            CTMLM.cloud_albedo(:lwp),
        ])
    end

    append!(eqs, extra_eqs)

    ds = processes_to_coupledodes(eqs, CTMLM;
        warn_default = false, name = :CTMLM, diffeq,
    )

    # The following parameters have been tuned to match observations
    # to an extend that was possible and reasonable.
    # The tuning is only really valid for a particular model setup
    # (the default model setup).

    if entrain == :Stevens2006
        set_parameter!(ds, :e_e, 0.8)
    elseif entrain == :Gesso2014
        set_parameter!(ds, :e_e, 0.5)
    end
    set_parameter!(ds, :d_c, 0.000875)
    set_parameter!(ds, :α_a, 0.275)
    try_set_parameter!(ds, :q_x_rate, 2.0)
    try_set_parameter!(ds, :s_x_0, 0.3)
    try_set_parameter!(ds, :T_FTR_0, 288.0) # this is the tropical TFTR in the paper
    try_set_parameter!(ds, :α_C_max, 0.6)

    for (p, val) in starting_parameters
        set_parameter!(ds, p, val)
    end

    return ds, eqs
end

sctebm_setup(input; kw...) = sctebm_setup(; kw..., input...)