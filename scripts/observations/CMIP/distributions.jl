using DrWatson
@quickactivate
using ClimateBase
using Statistics
include(srcdir("theme.jl"))
include(srcdir("cloud_boxes.jl"))
region = "Namibia"

# Data files
cmipdir(args...) = datadir("CMIP6", "direct_download", args...)

controlcodes = Dict(
    "ICON-ESM-LR" => "r1i1p1f1",
    "HadGEM3-GC31-LL" => "r4i1p1f3",
)
changecodes = Dict(
        "ICON-ESM-LR" => "r1i1p1f1",
        "HadGEM3-GC31-LL" => "r1i1p1f3",
)
yearcodes = Dict(
        "ICON-ESM-LR" => "201001-201412",
        "HadGEM3-GC31-LL" => "195001-201412",
)
yearchanges = Dict(
        "ICON-ESM-LR" => "201001-201912",
        "HadGEM3-GC31-LL" => "195001-199912",
)

vars = ["wap", "sfcWind"]

# load data while keeping region of interest and last year for all
box_corners = low_cloud_boxes_corners[region]

function specific_selection(var, file; remove_pressure = true, coorddim = false)
    if coorddim # icon files are in `CoordinateSpace` for example
        X = ncread(cmipdir(var*file), var; grid = CoordinateSpace())
        # wrap longitudes
    else
        X = ncread(cmipdir(var*file), var)
    end
    if remove_pressure && hasdim(X, Dim{:plev})
        # TODO: Make this into `Pre` dimension;
        X = X[Dim{:plev}(At(85000.0))]
    end
    # select last year of simulation
    t = size(X, Time)
    X = X[Time(t-11:t)]
    # wrap longitudes. Very important! Because
    # I've set up my boxes to have negative longitudes but
    # some CMIP data are set up to have only positive longitudes!
    # That is why we circshift the loaded data to have [-180, 180) longitudes!
    X = longitude_circshift(X)
    # finally select the appropriate cloud box
    selector = (
        Lat(Between(box_corners[2], box_corners[2]+low_cloud_box_lat)),
        Lon(Between(box_corners[1], box_corners[1]+low_cloud_box_lon)),
    )
    if coorddim
        selector = (Coord(selector...),)
    end
    return X[selector...]
end

fig = Figure()
l1 = l2 = l3 = nothing

for (jj, modelname) in enumerate(keys(controlcodes))

changecode = changecodes[modelname]
controlcode = controlcodes[modelname]
yearcode = yearcodes[modelname]
yearchange = yearchanges[modelname]
axs = axesgrid!(fig[1, jj], 2, 2;
    titles = ["U, $(modelname)", "Ω, $(modelname)"], ylabels = jj == 1 ? ["density", "ratio density"] : nothing)
change_file = "_Amon_$(modelname)_1pctCO2_$(changecode)_gn_$(yearchange).nc"
control_file = "_Amon_$(modelname)_historical_$(controlcode)_gn_$(yearcode).nc"


for (i, varname) in enumerate(("sfcWind", "wap"))

    V_control = specific_selection(varname, control_file; coorddim = modelname == "ICON-ESM-LR")
    V_change = specific_selection(varname, change_file; coorddim = modelname == "ICON-ESM-LR")
    V_ratio = V_change ./ V_control

    # we must do this simplification because of the very low Ω
    # values that may change sign
    if varname == "wap"
        V_ratio = vec(V_ratio[findall(o -> 0 ≤ o ≤ 2, V_ratio)])
    end

    # wind
    l1 = density!(axs[1, i], vec(V_control); label = "control", color = (COLORS[1], 0.5),
        strokecolor = COLORS[1], strokewidth = 2
    )
    vlines!(axs[1, i], mean(V_control); color = COLORS[1])

    l2 = density!(axs[1, i], vec(V_change); label = "change", color = (COLORS[2], 0.5),
        strokecolor = COLORS[2], strokewidth = 2, linestyle = LINESTYLES[2],
    )
    vlines!(axs[1, i], mean(V_change); color = COLORS[2], linestyle = LINESTYLES[2])

    # for the ratio we only plot quantiles to avoid excessive ranges
    q1, q2 = quantile(V_ratio, [0.025, 0.975])

    V_ratio = filter(v -> q1 < v < q2, V_ratio)


    l3 = density!(axs[2, i], vec(V_ratio); label = "ratio", color = (COLORS[3], 0.5),
    strokecolor = COLORS[3], strokewidth = 2
    )
    vlines!(axs[2, i], mean(V_ratio); color = COLORS[3])

    if i == 2
        xlims!(axs[1, i], 0, nothing)
        xlims!(axs[2, i], 0, nothing)
    end
end

hideydecorations!.(axs; grid = true, label = false)
ylims!.(axs, 0, nothing)


end

Legend(fig[1, 3], [l1, l2, l3], ["historical", "1% CO₂", "ratio"]; tellheight=false)

display(fig)
wsave(papersdir("figures", "cmip6"), fig)
