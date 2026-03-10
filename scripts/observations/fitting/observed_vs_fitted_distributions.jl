using DrWatson
@quickactivate
using Statistics
using Distributions
include(srcdir("theme.jl"))

# Stage 1: load observational data as input parameters
# this ensures we use exactly the same combinations
include(scriptsdir("observations", "ERA5_SCT", "generate_observations.jl"))

# %%
# Stage 2: Load simulated data
cooling = :q_x
invfix = :difference
Ld = :three_layer # has massive impact on the dynamics
ΔF_s = :ctrc # has massive impact on the dynamics
cloud_rad = :lwp
entrain = :Stevens2006
ftrgrad = :weak
input = @dict(cooling, invfix, cloud_rad, Ld, ΔF_s, entrain, ftrgrad)

foldergroup = ["sims", "observed_params"]

fitted_data = wload(datadir(foldergroup..., savename(input, "jld2")))["observables_distributions"]

# Stage 3: plot all the stuff together!
names = [
  :C => "",
  :SST => "K",
  :q_b => "g/kg",
  :CTRC => "W/m²",
  :Ld => "W/m²",
  :Lnet => "W/m²",
  :LHF => "W/m²",
  :RCT => "W/m²",
  :SHF => "W/m²",
  :ASW => "W/m²",
]

fitted_data[:BLH] = fitted_data[:z_b]

fig, axs = axesgrid(2, 5; size = (figwidth, 2.1*figheight))

for (specifier, data) in enumerate((outputs, fitted_data))
    for (k, name) in enumerate(names)
        i, j = Tuple(CartesianIndices((5,2))[k])
        ax = axs[j,i]
        x = data[name[1]]
        if specifier == 1 # observations must be filtered
            v = skipnan(x[valid_idxs])
        else
            v = vec(x)
        end
        alpha = specifier == 1 ? 0.1 : 0.2
        linestyle = specifier == 1 ? :solid : :dash
        den = density!(ax, v;
            color = (COLORS[k], alpha), strokecolor = COLORS[k], strokewidth = 3, linestyle
        )

        # means and stuff
        mS = round(mean(v); sigdigits = 3)
        preface = specifier == 1 ? "obs: " : "model: "
        textbox!(ax, preface*string(mS); valign = :top,
            halign = specifier == 1 ? :left : :right,
        )

        if specifier == 2
            hideydecorations!(ax)
            ylims!(ax, 0, nothing)
            xlabel = string(name[1])
            if !isempty(name[2])
                xlabel *= " [$(name[2])]"
            end
            ax.xlabel = xlabel
            # xlims!(ax, density_limits(name))
        end
    end
end


# Stage 4: add overarching axis
Legend(fig[0,:],
    [
        PolyElement(color = (:black, 0.1), strokewidth = 2, strokecolor = :black,),
        PolyElement(color = (:black, 0.2), strokewidth = 2, strokecolor=:black, linestyle = :dash)
    ],
    ["observations", "model"];
    tellwidth = false, tellheight = true, nbanks = 2
)

colgap!(fig.layout, 25)
display(fig)
wsave(papersdir("figures", "fitting_distributions"), fig)
