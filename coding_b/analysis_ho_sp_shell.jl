using ArgParse, DelimitedFiles, LinearAlgebra, Printf, Statistics, CSV, DataFrames
include("harm_osc.jl")

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "sample"
            help = "Number of steps in the simulation"
            required = true
            arg_type = Int
        "blocksize"
            help = "the number of points in a block"
            required = true
            arg_type = Int
        "therm"
            help = "the number of points to discard as thermalization"
            required = true
            arg_type = Int
        "--path", "-p"
            help = "the path where files are stored"
            default = joinpath(["..", "simulations_b"])
            required = false
            arg_type = String
        "--name", "-n"
            help = "the name of data file"
            default = "data.txt"
            required = false
            arg_type = String
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_cmd()
    path = parsed_args["path"]
    blocksize = parsed_args["blocksize"]
    therm = parsed_args["therm"]
    sample = parse_args["sample"]
    dfname = parse_args["name"]
    startp = @sprintf "ho_sp_sample=%.1e" sample
    paths = filter(startswith(startp), readdir(path))

    for fname in paths
        bend = findfirst("beta=", fname)[end]
        Nstr = findfirst("Nt=", fname)[1]
        Nend = findfirst("Nt=", fname)[end]
        Nt = parse(Int, fname[Nend+1:end-4])
        β = parse(Float64, fname[bend+1:Nstr-1])
        η = β/Nt
    
        local w = open(joinpath([path, fname]), "r") do io
            readdlm(io, header = true)
        end
        
        spacing=Vector{Int}()
        for el in w[2][2:end]
            append!(spacing, parse(Int, el[4:end]))
        end
    
        unique!(spacing)
        cormat = w[1][therm+1:end,:]
        datajack = zeros(size(cormat,1)÷blocksize, size(cormat,2))
    
        for i in 1:size(cormat,2)
            datajack[:,i] = JackKnife(cormat[:,i], blocksize)
        end
    
        corrected = datajack
        corrected[:,3:4:end].-=(corrected[:,1].^2)
        corrected = corrected[:,2:end]
        gaps = log.(corrected./circshift(corrected, (0,-4)))
        errs = std(gaps, dims = 1, corrected = false).*sqrt(size(datajack,1)-1)
        gaps = mean(gaps, dims = 1)
        gaps = gaps[:,1:end-4]./η
        errs = errs[:,1:end-4]./η

        data = DataFrame(
            sp = spacing[1:end-1],
            x = gaps[1:4:end],
            xv = errs[1:4:end],
            x2 = gaps[2:4:end],
            x2v = errs[2:4:end],
            x3 = gaps[3:4:end],
            x3v = errs[3:4:end],
            x4 = gaps[4:4:end],
            x4v = errs[4:4:end],
        )
        if !isdir(joinpath([path, "data"]))
            mkpath(joinpath([path, "data"]))
        end
        f1 = @sprintf "data_sample%.1eNt%ibeta%.1f.csv" sample Nt β
        touch(joinpath([path, "data", f1]))
        CSV.write(joinpath([path, "data", f1]), data) 
    end
end

main()
