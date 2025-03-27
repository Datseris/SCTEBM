using DrWatson
@quickactivate
using DataFrames
using Query
include(srcdir("theme.jl"))

# First, collect all outputs into a single dataframe
foldergroup = ["random_param_noco2"]
folder = datadir("sims", foldergroup...)
df = collect_results(folder)

# all simulation options
sim_options = [:ftrgrad, :cooling, :invfix, :ΔF, :w_m]

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

# Dataframes
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