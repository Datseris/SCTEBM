using DrWatson
@quickactivate
using ClimateBase
using ComplexityMeasures # for encoding into histograms
using Statistics
using Dates

file = datadir("CTRC", "CTRC_Zheng_GRL_2021.nc")

ncdetails(file)

tvec = vec(ncread(file, "time"))
hvec = vec(ncread(file, "UTC")) # hour time mark
lons = vec(ncread(file, "lon"))
lats = vec(ncread(file, "lat"))
vars = vec(ncread(file, "var_names"))
# we need to convert the data from singular vector to 3D lon-lat-time,
# and we will do so using monthly means and monthly stds.
# We will also only keep the following fields:
keep = ["CTRC", "SW_heat", "LW_cool", "LW_heat", "LW_CRE", "SW_heat_dmean",
"CF", "LWP", "COD", "CTT"]
var_idxs = [findfirst(isequal(v), vars) for v in keep]

# We will make the data in the same format I download ERA5 data at,
# which means only at a particular region (e.g., west of Namibia).
# Be sure it is exactly like ERA5 data!!!
region = "California"
if region == "Namibia"
    lons_era5 = -28.0f0:1.0f0:12.0f0
    lats_era5 = -5.0f0:-1.0f0:-25.0f0
elseif region == "California"
    lons_era5 = -155.0f0:1.0f0:-115.0f0
    lats_era5 = 35.0f0:-1.0f0:15.0f0
end

# to collect the indices of the bins that contain the data,
# we use the fantastic ComplexityMeasures.jl, arguably the
# most intuitive way to "encode" vectors into histogram bins
expand(r) = (s = step(r); ((first(r) - s/2):s:(last(r)+s/2)))
lons_bin = expand(lons_era5)
lats_bin = reverse(expand(lats_era5)) # requires sorted arrays!!!
binning = FixedRectangularBinning((lons_bin, lats_bin), true)
encoding = RectangularBinEncoding(binning)

sampling = "monthly"
if sampling == "monthly"
    tvec_proc = [Date(2014, i, 15) for i in 1:12]
elseif sampling == "daily"
    tvec_proc = Date(2014,1,1):Day(1):Date(2014,12,31)
elseif sampling == "weekly"
    tvec_proc = Date(2014,1,1):Day(7):Date(2014,12,31)
end

# %%
# Before going into any data, we obtain the indices of data
# that are within the region of interest and within the months.
# We will create a `ClimArray` whose elements are empty integer vectors,
# and then fill it with indices.
indices = Array{Vector{Int}}(undef, length.((lons_era5, lats_era5, tvec_proc)))
for i in eachindex(indices); indices[i] = Int[]; end
indices = ClimArray(indices, (Lon(lons_era5), Lat(lats_era5), Tim(tvec_proc)))

# search function to find time index
if sampling == "monthly"
    f = month
elseif sampling == "daily"
    f = dayofyear
elseif sampling == "weekly"
    f = week
end

# Right, now we go through the months and push! into the indices
for timi in eachindex(tvec_proc)
    # timi, loni, lati are the indices of the dimensions!
    # here we abuse the fact the the month/day number is also the index
    # thank god for 1-based indexing!
    timespan_idxs = findall(t -> f(t) == timi, tvec)

    for i in timespan_idxs

        j = encode(encoding, (lons[i], lats[i]))
        j == -1 && continue
        # convert linear index to cartesian
        ci = encoding.ci[encoding.li[j]]
        loni, lati = ci[1], ci[2]
        # note; the indices have reverse-order latitude,
        # so we transform the latitude index to its reverse
        lati = length(lats_era5) + 1 - lati

        # Maaaaan using ComplexityMeasures.jl made this so freaking easy!!!!
        # This is the type of code I was struggling with before:
        # if !(minimum(lons_search) < lons[i] < maximum(lons_search))
        #     continue
        # elseif !(minimum(lats_search) < lats[i] < maximum(lats_search))
        #     continue
        # end

        # loni = searchsortedfirst(lons_search, lons[i]) - 1
        # lati = searchsortedfirst(lats_search, lats[i]; rev = lats_rev)

        # if lats_rev
        #     lati = searchsortedfirst(lats_search, lats[i]; rev = true)
        # else
        #     lati = searchsortedlast(lats_search, lats[i])
        # end

        push!(indices[Tim(timi), Lon(loni), Lat(lati)], i)
    end
end


# %%
# Very good, now we can create the actual fields!
fields_mean = Any[
    ClimArray(zeros(length.(dims(indices))), dims(indices); name = n)
    for n in keep
]
fields_std = Any[
    ClimArray(zeros(length.(dims(indices))), dims(indices); name = n*"_std")
    for n in keep
]

# iterate through the fields, obtain the selection we want, get mean and std
for (fi, vi) in enumerate(var_idxs)
    field_mean = fields_mean[fi]
    field_std = fields_std[fi]
    # load the whole field
    X = ncread(file, "arrvar", ([vi], 1:length(tvec)))
    # then go through the space-time slices
    for (j, idxs) in enumerate(indices)
        Xsel = X[idxs]
        xmean = mean(Xsel)
        field_mean[j] = xmean
        field_std[j] = std(Xsel; mean = xmean)
    end
end

# %% test plot
using CairoMakie
F = fields_mean[1]
fig, ax, hm = heatmap(lons_era5, lats_era5, gnv(F[Tim(9)]))
cb = Colorbar(fig[1,2], hm)
fig

# %% Save!
globalattr = Dict(
    "description" => "Processing of CTRC_Zheng_GRL_2021.nc",
    "region" => region,
    "time average" => sampling,
    "original data source" => "https://zenodo.org/records/5218655",
)

ncwrite(datadir("CTRC", "CTRC_ZhengGRL_2014_$(region)_$(sampling).nc"),
    [fields_mean..., fields_std...]; globalattr
)
