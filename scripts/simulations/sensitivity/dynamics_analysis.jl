using DrWatson
@quickactivate
using DataFrames
using Query
include(srcdir("theme.jl"))

# First, collect all outputs into a single dataframe
foldergroup = ["random_param_noco2"]
folder = datadir("sims", foldergroup...)
df = collect_results(folder)

# fix missing entries for options that were not specified
replace!(df[!, "cloud_rad"], missing => :const)

# and finally sort dataframe by the most important aspects
sort!(df, ["cloud_rad"]; by = x -> string(x)[2], rev = false)
sort!(df, ["ΔF_s"]; by = x -> string(x)[2], rev = true)

# we don't need rMI values for this
select!(df, Not("rMI"))

# Gotta rename all names ending with '%' because it makes working
# with the macro-based querying systems impossible
function remove_ending_percentage!(df)
    n = names(df)
    problematic = findall(endswith("%"), n)
    renames = [p => p[1:end-1] for p in n[problematic]]
    rename!(df, renames...)
end
remove_ending_percentage!(df)


# %%
# find all with limit cycle using dataframes meta
x = @from row in df begin
    @where row.limitcycle > 0
    @select row
    @collect DataFrame
end

using DataFramesMeta

# Simulations that result in a limit cycle
x = @chain df begin
    DataFramesMeta.@rsubset(:limitcycle > 1)
    DataFramesMeta.@select($(sim_options...))
end
display(x)
# It appears that `w_m` is irrelevant for the presence of a limit cycle!
# If `q_x` cooling is used, we always get a limit cycle.
# If `sst_x` cooling is used, we also need `invfix` to be `difference`.
# Amazing!
# Note to reader: these comments do not make sense for the paper,
# as they refer to a super early version of the model!

# %% Any sorts of subsequent analysis here
using DataFramesMeta
sim_options = Symbol.(aspects)

# Find 10 most stable simulations
x = @chain df begin
    DataFramesMeta.@orderby(:diverged)
    DataFramesMeta.@rsubset(:diverged < 5)
    DataFramesMeta.@select($(sim_options...))
end
display(x)

# Okay it becomes clear that the most stable model configuration is with the stevens
# entrainment and the difference inversion condition.