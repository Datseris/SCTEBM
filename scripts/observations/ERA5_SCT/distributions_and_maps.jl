# This script uses loaded observational data and visualizes them as
# distributions and as maps
include("generate_observations.jl")

K = max(length(inputs), length(outputs))

##########################################################################################
# Distributions
##########################################################################################
csplit = 0.0
cidxs1 = findall(≥(csplit), outputs[:C])
cidxs2 = findall(<(csplit), outputs[:C])

idxs1 = intersect(cidxs1, valid_idxs)
idxs2 = intersect(cidxs2, valid_idxs)

figin = Figure(size = (3figwidth/4, figheight*3))
figout = Figure(size = (1200, 1000))

# Add units
units = Dict(
  :U => " [m/s]",
  :d_c => " [10⁻³]",
  :T_FTR => " [K]",
  :D => " [10⁻⁶ 1/s]",
  :α_a => "",
  :α_c => "",
  :S => " [W/m²]",
  :EIS => " [K]",
  :RH₊ => "",
)
inputs[:d_c] = inputs[:d_c] .* 1000
counter = 0
for (data, fig) in zip((inputs, outputs), (figin, figout))
    counter += 1
    for (k, name) in enumerate(keys(data))
        i, j = Tuple(CartesianIndices((3,5))[k])
        x = data[name]
        v = skipnan(x[idxs1]) # only at clouds we want
        mS = round(mean(v); sigdigits = 3)
        ax, den = density(fig[i, j], v;
            color = (COLORS[k], 0.5), strokecolor = COLORS[k], strokewidth = 3
        )
        textbox!(ax, string(mS))
        hideydecorations!(ax)
        ylims!(ax, 0, nothing)
        if !isempty(idxs2)
            v = skipnan(x[idxs2]) # only at clouds we want
            mC = round(mean(v); sigdigits = 3)
            density!(ax, v; label = string(mC), color = (COLORS[k], 0.25), linestyle = :dash, strokecolor = COLORS[k], strokewidth = 3)
        end

        if counter == 1
            ax.xlabel = string(name)*units[name]
        else
            ax.xlabel = string(name)
        end

        xlims!(ax, density_limits(name))
    end
end

colgap!(figin.layout, 25)
colgap!(figout.layout, 25)

# Save for paper without title
wsave(papersdir("figures", "observations_inputs"), figin)
# wsave(papersdir("figures", "observations_ouputs"), figout)

# Save with title for my reference
intit = "ERA5 - $(region) - inputs - $(date)"
outit = "ERA5 - $(region) - outputs - $(date)"
figuretitle!(figin, intit)
figuretitle!(figout, outit)

display(figin)
sleep(0.1)
display(figout)

# wsave(plotsdir("observations", intit*".png"), figin)
# wsave(plotsdir("observations", outit*".png"), figout)


##########################################################################################
# %% Maps
##########################################################################################
# So this doesn't work because spatial interruption due to D condition
idxs = setdiff(findall(isreal, f(M)), intersect(didxs, oidxs))

figin = Figure(size = (1200, 1000))
figout = Figure(size = (1200, 1000))

for k in 1:K # max number of in or out fields fields
    i, j = Tuple(CartesianIndices((3,5))[k])
    for (data, fig) in zip((inputs, outputs), (figin, figout))
        k > length(data) && continue
        X = data[k]
        x = deepcopy(f(X))
        x[idxs] .= NaN
        ax, hm = heatmap(fig[i, j], x;
            colormap = :dense
        )
        ax.xlabel = String(keys(data)[k])
    end
end

intit = "ERA5 - $(region) - inputs - $(date)"
outit = "ERA5 - $(region) - outputs - $(date)"
figuretitle!(figin, intit)
figuretitle!(figout, outit)

display(figin)
sleep(0.1)
display(figout)


# %% atmospheric profiles
# fig, ax, hm = heatmap(lon, lat, gnv(timemean(D[Pre(At(850))])))
fig, ax, hm = heatmap(lon, lat, gnv(timemean(α_c)); colormap = cgrad(:default, 9, categorical = true))
cb = Colorbar(fig[2,1], hm; vertical = false)
fig

# %%

idx = (Lon(At(5.0)), Lat(At(-15.0)), Tim(1))
# lines(T[idx...], Z[idx...])
lines(gnv(dims(Z, Pre)), gnv(Z[idx...]))
