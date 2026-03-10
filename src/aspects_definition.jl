# Make sure the aspects ordering is like in the paper table
# and make sure all match the options!
const ASPECTS_OPTIONS = OrderedDict(
    "ftrgrad" => [:none, :weak, :strong],
    "co2" => [1, 2, 3],
    "ΔF_s" => [:ctrc, :three_layer, :Gesso2014],
    "Ld" => [:three_layer, :fixed],
    "entrain" => [:Stevens2006, :Gesso2014],
    "cloud_rad" => [:const, :lwp],
)

"""
    aspects_to_option_idx(row)

Given a `DataFrame` row, or in general a named tuple
mapping aspects to their options, return the integer
(1, 2, 3) corresponding to that option
"""
function aspects_to_option_idx(row)
    sorted_aspects = collect(ASPECTS_OPTIONS)
    idxs = zeros(Int, length(ASPECTS_OPTIONS))
    # Ensure there is an informative error if any of these is nothing
    for i in eachindex(idxs)
        aspect, options = sorted_aspects[i]
        idx = findfirst(isequal(row[aspect]), options)
        if isnothing(idx)
            @show row aspect options
            error("index nothing")
        end
        idxs[i] = idx
    end
    # Same but in a one liner:
    # [findfirst(isequal(row[aspect]), options) for (aspect, options) in pairs(aspects)]
    return idxs
end
