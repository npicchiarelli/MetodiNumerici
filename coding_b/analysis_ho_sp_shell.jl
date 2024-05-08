using ArgParse, DelimitedFiles, LinearAlgebra, Plots, Printf, Statistics
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
    startp = @sprintf "ho_th_sample=%.1e" sample
    paths = filter(startswith(startp), readdir(path))
    size_list = []
    eta_list = []
    x = []
    xv = []
    x2 = []
    x2v = []
    K = []
    Kv = []
    H = []
    Hv = []

    for (i,fname) in enumerate(paths)
        Nt = parse(Int, fname[end-5:end-4])
        β = parse(Float64, fname[end-11:end-9])
        η = β/Nt

        local w = open(joinpath([path, fname]), "r") do io
            readdlm(io, header = true)
        end

        x_j = JackKnife(w[1][therm:end,1], blocksize)
        x2_j = JackKnife(w[1][therm:end,2], blocksize)
        K_j = JackKnife(w[1][therm:end,3], blocksize)
        H_j = 0.5.*x2_j.-K_j .+ 1/(2*η)
        
        push!(size_list, Nt)
        push!(eta_list, η)
        push!(x, mean(x_j))
        push!(x2, mean(x2_j))
        push!(K, mean(K_j))
        push!(H, mean(H_j))

        push!(xv, std(x_j, corrected = false).*sqrt(length(x_j)-1))
        push!(x2v, std(x2_j, corrected = false).*sqrt(length(x_j)-1))
        push!(Kv, std(K_j, corrected = false).*sqrt(length(x_j)-1))
        push!(Hv, std(H_j, corrected = false).*sqrt(length(x_j)-1))
    end

    w = open(joinpath([path, dfname]), "w") do io
        writedlm(io, ["size_list" "eta_list" "x" "xv" "x2" "x2v" "K" "Kv" "H" "Hv"], " ")
        writedlm(io, [size_list eta_list x xv x2 x2v K Kv H Hv], " ")
    end
end

