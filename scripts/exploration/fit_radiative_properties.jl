# For fitting emissivity cloud albedo based on expressions of old papers
using DrWatson
@quickactivate
include(srcdir("theme.jl"))

fig, axs = axesgrid(3, 1; size = (600, 400), xlabels = "LWP (g/m²)", ylabels = ["τ_C", "α_C", "ε_C", ], sharex = true)

lwp = 10 .^ range(-1, 3; length = 101)
epsilon = @. 1 - exp(-0.158*lwp)
tau = @. lwp^1.7/(10 + lwp)
g = 0.95
alpha = @.tau/(2/(sqrt(3)*(1 - g)) + tau)

lines!(axs[1], lwp, tau)
lines!(axs[2], lwp, alpha)
lines!(axs[3], lwp, epsilon)
axs[1].xscale = log10
axs[1].yscale = log10
axs[2].xscale = log10
axs[3].xscale = log10
axs[1].title = "fits of ε, τ, α from Stephens1978 and Datseris2021"

# wsave(plotsdir("observations", "fits_of_cloud_rad"), fig)

fig