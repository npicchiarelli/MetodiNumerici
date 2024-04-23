using ArgParse, Random

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
    
    Ntstr = @sprintf "%.1e" Nt
    η = β/Nt
    # initializing...
    lattice = zeros(Nt, Float64)

    # files management
    if !isdir(path)
        mkpath(path)
    end
    fname = "ho_"*"$Nt"*"$sample"*"$β"
    fr = joinpath([path, fname])
    if !isfile(fr)
        touch(fr)
    end
    # writing header
    open(f1, "w") do infile
        writedlm(file, ["x" "x2" "K" "c1" "c2" "c3" "c3c"], " ")
    end
    
end