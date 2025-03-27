# Provide a central function `ctmlm_setup` that creates
# a dynamical system representing the SCTEBM. Create any variant by utilizing
# keywords that map aspect to its options.

using ConceptualClimateModels
import ConceptualClimateModels.CloudToppedMixedLayerModel as CTMLM

"""
    ctmlm_setup(; aspect = option, ...) -> ds, eqs

Main function associated with the research article (Datseris 2025).
Given an input of SCTEBM aspect options, create and return the resulting dynamical system
and equations (that were given to `processes_to_coupledodes`).
The aspects are named slightly different from the Table 1 of the paper,
so it is best to inspect the source code of the function for their names.
"""
function ctmlm_setup(;
        dt = 0.05, diffeq = (adaptive = false, dt, verbose = false),
        extra_eqs = [], starting_parameters::Dict = Dict(),
        # The keywords are the aspects and their values are the options:
        invfix = :difference, # :difference or :temperature
        ftrgrad = :weak, # :none, :weak:, :strong, or number
        Ld = :three_layer, # :three_layer, :fixed, or other options, see function
        ΔF = :crc, # :ctrc, :three_layer, :Gesso2014
        co2::Int = 1, # 1 2 or 3 for increasing CO2 impact
        entrain = :Gesso14, # :Gesso2014 or :Stevens2006
        cooling = :q_x,
        w_m = false,
        cdversion = :sigmoid,
        ε_c = :fraction, T_c_version = :top,
        invdec = true, # enable decoupling augmentations with i𝒹
        fcrc = true, # fractional cloud radiative cooling, ∝ C
    )

    # Core equations
    eqs = [
        # Cloud equations
        CTMLM.cf_dynamic(cdversion),
        CTMLM.decoupling_variable(),
        CTMLM.cloud_layer_thickness(:Bolton1980),
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
        CTMLM.cloud_emission_temperature(T_c_version),
        CTMLM.bbl_emission_temperature(),

        # Radiation
        CTMLM.cloud_longwave_cooling(cloud_fraction = (ε_c ≠ :fraction && fcrc)),
        CTMLM.cloud_emissivity(ε_c),
        CTMLM.cloud_shortwave_warming(:insolation; cloud_fraction = fcrc),
        CTMLM.mlm_radiative_cooling(ΔF),
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
            CTMLM.𝒹_q ~ clamp((CTMLM.z_b*CTMLM.CLT/2750)^1.3, 0, 0.5),
            CTMLM.𝒹_s ~ 0.5*CTMLM.𝒹_q,
            CTMLM.Δ₊q ~ (1 - CTMLM.i_𝒟*CTMLM.𝒹_q)*(CTMLM.q₊ - CTMLM.q_b),
            CTMLM.Δ₊s ~ (1 - CTMLM.i_𝒟*CTMLM.𝒹_s)*(CTMLM.s₊ - CTMLM.s_b),
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
        set_parameter!(ds, :e_e, 0.7)
    elseif entrain == :Gesso2014
        set_parameter!(ds, :e_e, 0.5)
    end
    set_parameter!(ds, :c_d, 0.000875)
    set_parameter!(ds, :α_a, 0.275)
    # Actually tunable
    set_parameter!(ds, :q_x_rate, 2.0)
    set_parameter!(ds, :s_x_0, 0.2)
    set_parameter!(ds, :T_FTR_0, 288.0) # this is the tropical TFTR in the paper
    set_parameter!(ds, :α_C, 0.4)

    for (p, val) in starting_parameters
        set_parameter!(ds, p, val)
    end

    return ds, eqs
end

ctmlm_setup(input; kw...) = ctmlm_setup(; kw..., input...)