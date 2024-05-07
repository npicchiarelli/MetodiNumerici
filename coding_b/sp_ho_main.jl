using ArgParse, Dates, DelimitedFiles, Random, Printf
include("harm_osc.jl")

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "Nt"
            help = "the number of time steps"
            required = true
            arg_type = Int
        "sample"
            help = "Number of steps in the simulation"
            required = true
            arg_type = Int
        "beta"
            help = "beta = hbar*omega/(k_B T)"
            required = true
            arg_type = Float64
        "--path", "-p"
            help = "the path where files are stored"
            default = joinpath(["..", "simulations_b"])
            required = false
            arg_type = String
        "--metropolis", "-m"
            help = "if true, metropolis algo is used instead of heathbath"
            action = :store_true
        "--verbose", "-v"
            help = "if given, percentage of execution is printed"
            action = :store_true
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_cmd()
    Nt = parsed_args["Nt"]
    sample = parsed_args["sample"]
    β = parsed_args["beta"]
    path = parsed_args["path"]
    metropolis = parsed_args["metropolis"]
    verbose = parsed_args["verbose"]
    
    sstr = @sprintf "%.1e" sample
    η = β/Nt
    # initializing...
    lattice = zeros(Float64, Nt)
    acc = 0.

    # simulation parameters
    orsteps = 5
    Δ = 10.0*η
    measevery = 1000

    # files management
    if !isdir(path)
        mkpath(path)
    end
    fname = "ho_sp"*"Nt=$Nt"*"sample=$sstr"*"beta=$β"*".txt"
    fr = joinpath([path, fname])
    if !isfile(fr)
        touch(fr)
    end
    # writing header
    cs = permutedims(["c$(Int(di))_$ci" for ci in 1:(Nt÷4)+1 for di in 1:4])
    open(fr, "w") do infile
        writedlm(infile, cs, " ")
    end
    
    start = now()
    datafile = open(fr, "a")
    for iter in 1:sample
        for r in LinearIndices(lattice)
            if rand() < .5
                if metropolis
                    acc+=metropolis!(lattice, r, Δ, η)
                else
                    acc+=heathbath!(lattice, r, η)
                end
            else
                for _ in 1:orsteps
                    acc+=overrelax!(lattice, r, η)
                end
            end
        end

        if iter%measevery == 0
            corr = []
            for δt in 0:Nt÷4
                append!(corr, correlators(lattice, δt))
            end
            writedlm(datafile, [permutedims(corr)], " ")
        end
        if verbose && iter % (sample÷100) == 0
            print("$((100*iter÷sample))% \r")
        end
    end
    close(datafile)
    elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
    println("\n$(round(now(), Dates.Second));\nNₜ = $Nt, elapsed time $(elapsed)\n")
end
main()