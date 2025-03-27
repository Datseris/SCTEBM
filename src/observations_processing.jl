using ClimateBase

function reprocess_to_daily(datapath, nameprefix, D, year)

    filegen() = "$(D)D_$(nameprefix)_$(sampling)_$(year)"
    filegen(month) = "$(D)D_$(nameprefix)_$(sampling)_$(year)_$(month)"

    # step 1: process each month to daily means
    for month in 1:12
        filename = filegen(month)
        file = joinpath(datapath, filename)*".nc"

        variables = nckeys(file)
        filter!(v -> v ∉ ("longitude", "latitude", "level", "time"), variables)

        processed_vars = Any[]
        for v in variables
            X = ncread(file, v)
            Xd = dailyagg(X)
            push!(processed_vars, Xd)
        end

        filename = filegen(month)*"_proc"
        ncwrite(joinpath(datapath, filename)*".nc", processed_vars)
    end

    # step 2: combine all months into one file
    files = filegen.(1:12) .* "_proc"
    files = [joinpath(datapath, file)*".nc" for file in files]
    ncd = NCDataset(files; aggdim = "time")
    file = joinpath(datapath, filegen()*".nc")
    NCDatasets.write(file, ncd)

    # step 3: delete all other files we used to process!!!
    for month in 1:12
        filename = joinpath(datapath, filegen(month))
        rm(filename*".nc")
        rm(filename*"_proc.nc")
    end
end
