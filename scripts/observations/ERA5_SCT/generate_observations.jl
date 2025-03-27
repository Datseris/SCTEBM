# Creates the two dictionaries `inputs, outputs`
# that is used in other scripts
# that utilize observational data. It always creates
# a cloud fraction field `C` and an ocean mask field
using DrWatson
@quickactivate
using ClimateBase
using Statistics
using CairoMakie; CairoMakie.activate!()
include(srcdir("theme.jl"))
include(srcdir("observations_analysis.jl"))

year = 2014
region = "Namibia"
sampling = "monthly"
input = @strdict(year, region, sampling)
fname(N) = "$(N)D_$(region)_$(sampling)_$(year).nc"
file2 = datadir("ERA5", fname(2))
file3 = datadir("ERA5", fname(3))

# we always need clouds, land mask, height, temperature, and RH
C = ncread(file2, "lcc") # low cloud fraction only
lon, lat = gnv.(dims(C, (Lon, Lat)))
M = ncread(file2, "lsm")
skipnan(v) = filter(!isnan, collect(skipmissing(v)))
Z = ncread(file3, "z") ./ 9.81
T = ncread(file3, "t")
Tsfc = ncread(file2, "t2m")
SST = ncread(file2, "sst")

# Then, we load all other fields we will need
LHF = -ncread(file2, "mslhf")
SHF = -ncread(file2, "msshf")
BLH = ncread(file2, "blh")
RH = ncread(file3, "r") ./ 100
D = ncread(file3, "d") ./ 1e-6
U = ncread(file2, "si10")

# This function computes the temperature inversion directly by finding
# the heights at which it increases with height. But the results are
# really crap, it is almost always 0. Therefore it is ignored in favour of EIS
# IS, IE, TS, TE = temperature_inversion(T, Z)
# TI = TE .- TS

# We need to obtain a representative relative humidity
RH_b = 0.8 # this is what Wood2006 uses
# and the average of this is almost identical:
RH_b = dropagg(mean, RH[Pre(Between(1000, 950))], Pre)

EIS = estimated_inversion_strength(Tsfc, T[Pre(At(700))], T[Pre(At(850))]; RH = RH_b)
LCL = lifting_condensation_level(Tsfc, RH_b)
CBH = ncread(file2, "cbh")

CLT = @. clamp((BLH - CBH)/(BLH), 0, 1)
# CLT = @. clamp((BLH - LCL)/(BLH), 0, 1)

# We modify this so that it also has liquid water
Q = (ncread(file3, "q") .+ ncread(file3, "clwc")) .* 1000
q_b = dropagg(mean, Q[Pre(Between(1000, 950))], Pre)

# From LHF, moisture and surface temperature, I can estimate the c_d coefficient
# q_saturation from Stevens 2006.
function q_saturation(T)
    ℓ_v = 2.53e6 # latent heat of vaporization (J/gr) # notice the gram units! coz we use `q` in gr/kg!
    p₀ = 101780.0 # reference pressure at the sea surface in Pa
    Rv = 461.0 # gas constant water vapor (J/K/kg)
    Rd = 287.0 # gas constant dry air (J/K/kg)
    e0 = 610.78 # reference water vapor partial pressure at 273.15 K in Pa = J/m^3
    psat = e0 * exp(-ℓ_v/Rv * (1/T - 1/273.16))
    return Rd/Rv * psat / (p₀ - psat)*1e3 # multiply with 1e3 due to my chosen units!
end
q₀ = q_saturation.(SST)

c_d = LHF ./ (q₀ - q_b) / (2.53e3) / 1.2 ./ U # constant density

# longwave
Ld = ncread(file2, "msdwlwrf")
Lnet = -ncread(file2, "msnlwrf") # we make all quantities positive
L₀ = Ld + Lnet
# ε_sfc = @. L₀/Tsfc^4/5.6703744e−8 # this is almost perfectly 1, but actually a bit > 1.

# albedos
S = ncread(file2, "mtdwswrf") # insolation
F_s_⬇ = ncread(file2, "msdwswrf")
F_s_⬇_K = ncread(file2, "msdwswrfcs")
F_s_net = ncread(file2, "msnswrf")
F_s_⬆ = F_s_⬇ .- F_s_net
α_s = F_s_⬆ ./ F_s_⬇
α_t = 1 .- (F_s_⬇ ./ S) # total, atmosphere + clouds
α_a = 1 .- (F_s_⬇_K ./ S)
α_c = @. 1 - (1 - α_t)/(1 - α_a)
α_c = α_c ./ C # normalize by cloud fraction

# select height for above BL
# RH₊ = RH[Pre(At(700))]
RH₊ = dropagg(mean, RH[Pre(Between(700, 500))], Pre)

# Unfortunately the resulting distributions depend strongly
# on the height where we obtain D, due to the filtering condition
# of D > 0... The original plan was to get D at pressures 900 or greater,
# which should be above the BL. However the typical values of D there are
# incredibly low, peaking at 1-2, which is so much lower than the typical values
# used in MLM studies that range in 4-6.
# So I decided to go with the value of D at the surface which is well within
# the range used in MLM studies.
D = D[Pre(At(1000))]
# we then assume the inversion height to be at:
Z₊ = dropagg(mean, Z[Pre(At(900))], (Tim,))
T₊ = T[Pre(At(900))]
T_FTR = T[Pre(At(850))]

# Radiative cooling
filectrc = datadir("CTRC", "CTRC_ZhengGRL_$(year)_$(region)_$(sampling).nc")
CTRC = ncread(filectrc, "CTRC")
# We need to scale CTRC by the cloud fraction to have all sky CTRC (discussions with Zheng)
CFZheng = ncread(filectrc, "CF")
CTRC = - CTRC .* CFZheng # - because we want it positive definite

ASW = F_s_net # Absorbed

# I am not plotting a_s, TI, z_+, as I do not use them anywhere
inputs = @dict D EIS RH₊ U T_FTR c_d α_a α_c S
outputs = OrderedDict(
    :C => C, :SST => SST, :q_b => q_b, :BLH => BLH,
    :CTRC => CTRC, :Ld => Ld, :Lnet => Lnet, :LHF => LHF,
    :CLT => CLT, :SHF=> SHF, :T₊ => T₊, :ASW => ASW,
)

# Finally we decide how to reduce the time dimension. For example:
# f(X) = timemean(X) # to include all time
# or
ti = 9 # month index
date = dims(C, Tim)[ti]
f(X) = gnv(X[Tim(ti)]) # to select a particular month

# and remap all fields to this
M = f(M)
inputs = OrderedDict(k => f(v) for (k, v) in inputs)
outputs = OrderedDict(k => f(v) for (k, v) in outputs)

# Lastly, collect all the indices that are valid given our approach
oidxs = findall(m -> m < 0.2, M) # ocean indices
# and also select points where there is subsidence
didxs = findall(>(0), inputs[:D])
# and also ignore points with 0 CLT
cltidxs = findall(clt -> !ismissing(clt) && clt != 0, outputs[:CLT])
# and combine everything
invalid_idxs = intersect(didxs, oidxs, cltidxs)