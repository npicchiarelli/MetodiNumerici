using ArgParse, CSV, DataFrames, Dates, Printf
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
        "beta"
            help = "the reverse temperature"
            required = true
            arg_type = Float64
        "metropolis"
            help = "if true, metropolis algo is used instead of heathbath"
            default = false
            arg_type = Bool
    end
    return parse_args(s)
end

function main()
    Q = 4
    Δ = 1
    parsed_args = parse_cmd()
    path = parsed_args["path"]
    Nt = parsed_args["Nt"]
    Ntstr = @sprintf "%.1e" Nt
    L = parsed_args["L"]
    β = parsed_args["beta"]
    metropolis = parsed_args["metropolis"]

    if isdir(path) == false
        mkpath(path)
    end

    Nupdates = Nt*L*L
    println("You are simulating a clock model with L=$L, Nt=$Nt, total number of updates = $Nupdates.\nThe files will be saved in $path")

    pdict = init_prob_dict(Q, β)
    E = Vector{Float64}(undef, Nt)
    m = Vector{ComplexF64}(undef, Nt)
    
    lattice = zeros(Int, (L,L)) 
    fname = "clock_"*"Nt="*"$Ntstr"*"L=$L"*"beta=$β"*".csv"
    f1 = joinpath([path, fname])
    
    touch(f1)
    acc = 0
    start = now()
    for nt in 1:Nt
        for idx in CartesianIndices(lattice)
            if metropolis
                rx = idx[1]
                ry = idx[2]
                acc += metropolis!(lattice,rx,ry,Δ,L,β) 
            else
                acc += heathbath!(lattice, idx, pdict, L)
            end
        end
    E[nt] = energy(lattice, Q, L)
    m[nt] = magnetization(lattice, Q, L)
    if nt % (Nt÷100) == 0
        print("$((100*nt/Nt))%...")
    end
    end
    
    datafile = open(f1)
    data = DataFrame(
        E = E,
        m = m
    )
    CSV.write(f1, data)
    close(datafile)
    local elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
    println("\n$(round(now(), Dates.Second));\nβ = $β, elapsed time $(elapsed)\n")
end

main()