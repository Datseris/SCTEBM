# Script for playing around with different equations and parameters for CTMLM
# There are two ways to specify the model equations;
# 1. use the `ctmlm_setup` function (this file) that makes same setups as in the paper
# 2. directly create your own equations (utilizing the default equations in defaults.jl),
# which is in the file `playground_equations.jl`

# From there, you specify what observables to visualize, and for which
# parameters to create sliders for. Both of these options are given
# after model initialization, which allows one to change change parameter
# sliders or observables visualized without restarting the model!

using DrWatson
@quickactivate
using DynamicalSystems
include(srcdir("ctmlm_setups.jl"))
include("playground_helpers.jl")

@parameters p = 5.0 e_m = 0.3 SST_X_0 = 0.0 q_x_rate = 2 β₋ = 0.3
@variables fake(t) = 0.5
@variables obs(t)
@variables clt2(t)
@variables w_e2(t)

eqs = [
    # clouds
    CTMLM.cf_dynamic(:sigmoid),
    CTMLM.decoupling_parameter(),
    # TimeDerivative(CTMLM.C, 0), # use this to make cloud fraction a fixed number
    CTMLM.CLT ~ CTMLM.cloud_layer_thickness(:Bolton1980),
    CTMLM.cloud_emissivity(:fraction),
    # Boundary layer standard stuff
    CTMLM.sst_dynamic(),
    CTMLM.mlm_dynamic(),
    CTMLM.mlm_q₊(:relative),
    CTMLM.mlm_s₊(:difference;
        cloud_effect = false, CO2_effect = false,
    ),
    # Temperatures
    CTMLM.free_troposphere_emission_temperature(:weak; add_co2 = false),
    CTMLM.cloud_emission_temperature(:top),
    CTMLM.bbl_emission_temperature(),
    # Radiation
    CTMLM.cloud_longwave_cooling(
        # cloud emissivity decides this because otherwise we don't have cloud
        # fraction affecting the longwave component
        cloud_fraction = false
    ),
    CTMLM.cloud_shortwave_warming(:insolation),
    CTMLM.bbl_radiative_cooling(:ctrc),
    CTMLM.downwards_longwave_radiation(:three_layer),
    CTMLM.entrainment_velocity(:Stevens2006; w_e = CTMLM.w_e, use_augmentation = true),

    # Cooling
    CTMLM.q_x ~ q_x_rate*max(CTMLM.q_saturation(CTMLM.SST)/CTMLM.q_saturation(290) - 1, 0),
    # CTMLM.q_x ~ p*max(CTMLM.q_saturation(CTMLM.SST)/CTMLM.q_saturation(290) - 1, 0),
    # CTMLM.q_x ~ p*max(CTMLM.q_b - 15, 0)^3,
    # CTMLM.q_x ~ p*max(CTMLM.q_saturation(CTMLM.s_b)/CTMLM.q_saturation(300) - 1, 0),
    ParameterProcess(CTMLM.s_x, 0.0),

    # CTMLM.ζ ~ CTMLM.CLT*β₋*(1 - CTMLM.C)/0.7/2,
    CTMLM.ζ ~ CTMLM.i_𝒟*CTMLM.CLT*β₋/2,
    # CTMLM.bbl_Δ₊sᵥ(),
    # For demonstration:
    # TimeDerivative(CTMLM.ζ, CTMLM.CLT*β₋*(1 - CTMLM.C)/2, 0.25),
    # CTMLM.w_e ~ 1e-3(fake*(CTMLM.SST/300) + 0.5),
    # fake ~ cos(t)^2,

    # CTMLM.LWP ~ CTMLM.liquid_water_path(CTMLM.T_t, CTMLM.CLT, CTMLM.z_b, CTMLM.s_b, CTMLM.q_b),

    # check additional variables:
    # clt2 ~ CTMLM.cloud_layer_thickness(:Bolton1980),
    # CTMLM.entrainment_velocity(:LL96DG14; w_e = w_e2, e_e = 0.7),

    # CTMLM.entrainment_velocity(:DalGesso2014; w_e = w_e2, use_augmentation = true),

    # Here I am testing enabling the decoupling
    CTMLM.𝒹_q ~ clamp((CTMLM.z_b*CTMLM.CLT/2750)^1.3, 0, 0.5),
    CTMLM.𝒹_s ~ 0.5*CTMLM.𝒹_q,
    # You can re-form analytically the expressions
    # from the `𝒹` definitions to my existing inversion as:
    CTMLM.Δ₊s ~ (1 - CTMLM.i_𝒟*CTMLM.𝒹_s)*(CTMLM.s₊ - CTMLM.s_b),
    CTMLM.Δ₊q ~ (1 - CTMLM.i_𝒟*CTMLM.𝒹_q)*(CTMLM.q₊ - CTMLM.q_b),
    # But still even with those the Sc state has larger height than the Cu.
    # Only if one changes the ΔF to not be proportional to the cloud fraction
    # one restores the "correct" balance of heights... E.g.:
    # CTMLM.ΔF ~ 40.0,
    # fake ~ CTMLM.Cᵢ,
    CTMLM.CRC ~ CTMLM.CRClw - CTMLM.CRCsw,

    # obs ~ CTMLM.cloud_layer_thickness(:Bolton1980)
]

ds = processes_to_coupledodes(eqs, CTMLM;
    warn_default = false, name = :CTMLM,
    diffeq = (adaptive = false, dt = 0.05)
)

set_parameter!(ds, :q_x_rate, 2)

# %% Configure and launch GUI
# the `GUI_obs` can be a premade configuration or a vector of things to observe
# which can include arbitrary symbolic expressions or functions as per DynamicalSystems.jl
# GUI_obs = :temperature
GUI_obs = [:LHF, :CLT, :ΔF, :CTRC, CTMLM.CRClw/CTMLM.C]
# The `GUI_par` can only be a vector of symbols
GUI_par = [:U, :D, :δ_Δ₊T, :α_C]

fig, dsobs = ctmlm_gui(ds, GUI_par, GUI_obs)
