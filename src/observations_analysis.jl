# File including functions used to analyze observation data
using ClimateBase
using OrderedCollections

module Constants
const Lhvap = 2.5e6    # Latent heat of vaporization (J / kg)
const Lhsub = 2.834e6   # Latent heat of sublimation (J / kg)
const Lhfus = Lhsub - Lhvap  # Latent heat of fusion (J / kg)
const cp = 1004.0     # specific heat at constant pressure for dry air (J / kg / K)
const Rd = 287.0         # gas constant for dry air (J / kg / K)
const R_over_cp = Rd / cp
const Rv = 461.5       # gas constant for water vapor (J / kg / K)
const cpv = 1875.0     # specific heat at constant pressure for water vapor (J / kg / K)
const eps = Rd / Rv
const g = 9.8          # gravitational acceleration (m / s**2)
const rho_w = 1000.0    # density of water (kg / m**3)
const cw = 4181.3      # specific heat of liquid water (J / kg / K)
end

###########################################################################################
# Albedo
###########################################################################################
"""
    effective_cloud_albedo(EBAF)

Estimate the _energetically consistent effective cloud albedo_ as in
Datseris & Stevens 2021 by loading all necessary fields from the given .nc file.
Assumes there is a single .nc file with all CERES EBAF fields.
"""
function effective_cloud_albedo(EBAF = datadir("observations", "CERES_EBAF_gea250.nc"))
    R = ncread(EBAF, "toa_sw_all_mon")
    K = ncread(EBAF, "toa_sw_clr_t_mon")
    I = ncread(EBAF, "solar_mon")
    F = ncread(EBAF, "cldarea_total_daynight_mon") ./ 100
    T = ncread(EBAF, "cldtau_total_day_mon")
    # Correct optical depth
    T = ClimateBase.sinusoidal_continuation(T, [1, 2.0]; Tmin = 0)

    # use `\:arrow_down:` and tab for arrows
    F_s_⬆ = ncread(EBAF, "sfc_sw_up_all_mon")
    F_s_⬇ = ncread(EBAF, "sfc_sw_down_all_mon")
    F_s_⬆_K = ncread(EBAF, "sfc_sw_up_clr_t_mon")
    F_s_⬇_K = ncread(EBAF, "sfc_sw_down_clr_t_mon");

    l = size(F_s_⬆, Time)
    argsall = timemean.((I[Time(1:l)], R[Time(1:l)], F_s_⬆, F_s_⬇))
    argsclr = timemean.((I[Time(1:l)], K[Time(1:l)], F_s_⬆_K, F_s_⬇_K))
    α_ATM, α_SFC = surface_atmosphere_contributions(argsall...) # all sky
    α_ATM_K, α_K_SFC = surface_atmosphere_contributions(argsclr...) # clear sky

    # This is the normalization done to have energetically consistent albedo value
    # according to Donohoe & Battisti model
    C2_timeavg = α_ATM .- α_ATM_K
    a = C2_timeavg; τ = timemean(T); f = timemean(F)
    x = -2a ./ (τ .* a - f .* τ)
    g_normalized = clamp.(1 .- x ./ eltype(τ)(√3), 0, 1)

    C = lacis_formula_cloud_albedo(F, T, g_normalized)
    return ClimArray(C; name = "α_cloud", attrib = Dict("standard_name" => "cloud albedo (CERES)"))
end


"""
    surface_atmosphere_contributions(I, F_toa_⬆, F_s_⬆, F_s_⬇) → α_ATM, α_SFC

Calculate the atmospheric and surface _contributions_ of the planetary albedo, so that
the TOA albedo is `α = α_ATM + α_SFC`, using the
simple 1-layer radiative transfer model by Donohoe & Battisti (2011) or G. Stephens (2015).
Stephens' formulas are incorrect and I have corrected them!
"""
function surface_atmosphere_contributions(I, F_toa_⬆, F_s_⬆, F_s_⬇)
    R = F_toa_⬆ ./ I   # planetary albedo (system reflectance)
    T = F_s_⬇ ./ I     # system transmisttance
    α = F_s_⬆ ./ F_s_⬇ # surface albedo

    # Formulas by Stephens, which are wrong!
    # τ = @. (1 - α*R) / (1 - (α^2) * (T^2)) # atmosphere transmittance
    # r = @. R - τ*α*T # atmospheric contribution to albedo

    # My calculations:
    r = @. (R - α*T^2) / (1 - (α*T)^2) # atmospheric contribution to albedo
    t = @. T*(1 - r*α)                 # atmospheric transmittance
    s = @. (α*t^2) / (1 - r*α)         # surface contribution to albedo
    return r, s
end

# This implements the formula from Lacis, multiplied with cloud fraction
function lacis_formula_cloud_albedo(f, τ, g = 0.9; frequencies = [1.0, 2.0])
    p = eltype(τ).(@. (sqrt(3))*(1 - g))
    R = p .* τ ./ (2 .+ p .* τ)
    cloud_frac = any(x -> x > 1, f) ? f ./ 100 : f # ensure in (0, 1)
    eca = cloud_frac .* R
    gstr = g isa Real ? string(g) : "array"
    return ClimArray(eca.data, eca.dims; name = "effective cloud albedo, g = $(gstr)")
end

###########################################################################################
# Inversion
###########################################################################################
"""
    estimated_inversion_strength(Tsfc, T700 [, T850]; Psfc = 1017.8, RH = 0.8) → EIS
Compute the Estimated Inversion Strength or EIS, from eq.(4) of [^Wood2006],
given `Tsfc` the surface temperature and `T700` the air temperature at 700 hPa
(in Kelvin).

EIS is a normalized measure of lower tropospheric stability acccounting for
temperature-dependence of the moist adiabat.

[^Wood2006]: On the relationship between stratiform low cloud cover and lower-tropospheric stability,
             Wood and Bretheron, [DOI: 10.1175/JCLI3988.1](https://journals.ametsoc.org/view/journals/clim/19/24/jcli3988.1.xml)
"""
function estimated_inversion_strength(Tsfc, T700, T850 = (Tsfc+T700)/2; Psfc = 1017.8, RH = 0.8)
    # Lower Tropospheric Stability (θ700 - θ0)
    LTS = lower_tropospheric_stability(Tsfc, T700, Psfc)
    # Default 80% relative humidity to compute LCL is appropriate for marine boundary layer
    LCL = lifting_condensation_level(Tsfc, RH)
    #  Γm = -dθ/dz is the rate of potential temperature decrease along the moist adiabat in K / m
    Γm = @. (Constants.g/Constants.cp*(1.0 - (1.0 + Constants.Lhvap*qsat(T850,850) / Constants.Rd / T850) /
         (1.0 + Constants.Lhvap^2 * qsat(T850, 850)/ Constants.cp/Constants.Rv/T850^2)))
    # Assume exponential decrease of pressure with scale height given by surface temperature
    z700 = @. (Constants.Rd*Tsfc/Constants.g)*log(Psfc/700.0)
    EIS = @. LTS - Γm*(z700 - LCL)
    return EIS
end

"""
    lower_tropospheric_stability(Tsfc, T700, Psfc = 1017.8)
By definition the difference of the potential temperature at 700hPa minus that at surface
pressure, by default = 1000.0 hPa.
"""
function lower_tropospheric_stability(Tsfc, T700, Psfc = 1017.8)
    @. potential_temperature(T700, 700) .- potential_temperature(Tsfc, Psfc)
end

function potential_temperature(T, P)
    @. T*(1017.8/P)^Constants.R_over_cp
end

"""
    lifting_condensation_level(T, RH) → LCL (meters)
Compute the Lifiting Condensation Level (LCL) for a given temperature and relative humidity.

This is height (relative to parcel height with temperature T) at which the parcel would become saturated
during adiabatic ascent. It is based on an approximate formula from Bolton
(1980 MWR) as given by Romps (2017 JAS).
For an exact formula see Romps (2017 JAS), doi:10.1175/JAS-D-17-0102.1
"""
function lifting_condensation_level(T, RH)
    Tadj = @. T - 55.0  # in Kelvin
    return @. Constants.cp/Constants.g*(Tadj - (1/Tadj - log(RH)/2840.0)^(-1))
end

"""
    qsat(T,P) → qₛ (dimensionless)
Compute saturation specific humidity as function of temperature (in Kelvin)
and pressure (in hPa).
"""
function qsat(T,p)
    es = clausius_clapeyron(T)
    q = @. Constants.eps * es / (p - (1 - Constants.eps) * es )
    return q
end

"""
    clausius_clapeyron(T) → es (hPa)
Compute saturation vapor pressure as function of temperature T (in Kelvin).

Formula from Rogers and Yau "A Short Course in Cloud Physics" (Pergammon Press), p. 16
claimed to be accurate to within 0.1% between -30degC and 35 degC
Based on the paper by Bolton (1980, Monthly Weather Review).
"""
function clausius_clapeyron(T)
    Tcel = T .- 273.15
    es = @. 6.112 * exp(17.67*Tcel/(Tcel+243.5))
    return es
end

"""
    temperature_inversion(T, z) → IS, IE, TS, TE

Return the height of the start and stop of temperature inversion as well as
the temperature inversion itself.
The inversion is detected as the first point where `T` first
increases with height and ends at the last point where `T` increases with height.
(assumed that given inputs are T(z), z).
Return `NaN, NaN, 0, 0` if there is no point where `T` increases with height.
"""
function temperature_inversion(T::AbstractVector, z::AbstractVector)
    length(T) ≠ length(z) && error("unequal length")
    si = sortperm(z)
    z = z[si]
    T = T[si]
    inv = false
    is = 1 # start and stop of inversion
    for i in 2:length(T)
        if inv == false
            if T[i] > T[i-1]
                is = i
                inv = true # we are in the inversion range
            end
        else
            if T[i] < T[i - 1]
                return (z[is], z[i-1], T[is], T[i-1])
            end
        end
    end
    return (NaN, NaN, 0, 0) # if we didn't find inversion
end
function temperature_inversion(T::ClimArray, Z::ClimArray)
    odims = otherdims(T, Pre())
    IS = ClimArray(zeros(length.(odims)), odims; attrib = Dict("long_name" => "Inversion start"))
    IE = ClimArray(zeros(length.(odims)), odims; attrib = Dict("long_name" => "Inversion end"))
    TS = ClimArray(zeros(length.(odims)), odims; attrib = Dict("long_name" => "Temperature inversion start"))
    TE = ClimArray(zeros(length.(odims)), odims; attrib = Dict("long_name" => "Temperature inversion end"))
    for j in collect(otheridxs(T, Pre())) # get height profiles
        t = gnv(T[j...])
        z = gnv(Z[j...])
        is, ie, ts, te = temperature_inversion(t, z)
        IS[j...] = is
        IE[j...] = ie
        TS[j...] = ts
        TE[j...] = te
    end
    return IS, IE, TS, TE
end
