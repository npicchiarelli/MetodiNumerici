using ArgParse, CSV, DataFrames, Dates, Printf, Statistics
include("clock.jl")


function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "path"
            help = "the path where files are stored"
            default = joinpath(["..", "simulations_a"])
            required = true
            arg_type = String
        "Nt"
            help = "the number of steps in the simulation"
            required = true
            arg_type = Int
        "L"
            help = "the lattice dimension"
            required = true
            arg_type = Int
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_cmd()
    path = parsed_args["path"]
    Nt = parsed_args["Nt"]
    Ntstr = @sprintf "%.1e" Nt
    L = parsed_args["L"]
    startp = "clock_"*"Nt="*"$Ntstr"*"L=$L"
    paths = filter(startswith(startp), readdir(path))

    f_w = path*"data"*"L=$L"*"Nt=$Ntstr"*".csv"
    touch(f_w)
    tosave = DataFrame([[],[],[],[],[],[],[],[],[],[],[],[],[],[]], ["beta", "m", "e","e_v","m_abs","m_abs_v","m_abs2","m_abs2_v","spec_heat","spec_heat_v","susc","susc_v","bind","bind_v"])

    for (i,fname) in enumerate(paths)
        ffile = open(joinpath([path, fname]))
        β = parse(Float64, fname[length(startp)+length("beta=")+1:end-4])

        df= CSV.read(ffile, DataFrame, types = [Float64, ComplexF64])

        energ = df[!,:E]
        magn = df[!,:m]

        therm = 1
        start = now()

        
        means, stds = JackKnife(energ[therm:end],magn[therm:end],100, L)
        inter = ([means' stds']'[:])'
        data = [β abs(mean(magn)) inter]

        push!(tosave, data)
        println("$L: β = $β,$i/$(size(paths,1))")
        close(ffile)
    end 
    CSV.write(f_w, tosave)
end

main()