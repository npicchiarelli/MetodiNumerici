using CSV, DataFrames, Dates
include("clock.jl")

Q = 4
mant = 5
expo = 5
Nt = Int(mant*10^expo)
Δ = 1
metropolis = false

path = joinpath(["..", "simulations_a"])

if isdir(path) == false
    mkpath(path)
end

size_list = [20 40 50]
for L in size_list

    Nupdates = Nt*L*L
    println("You are simulating a clock model with L=$L, Nt=$Nt, total number of updates = $Nupdates.\nThe files will be saved in $path")

    betarray = round.(LinRange(0.86,0.88,15), digits = 4)
    start_all = now()

    for (i, β) in enumerate(betarray[1:end])
        pdict = init_prob_dict(Q, β)
        E = Vector{Float64}(undef, Nt)
        m = Vector{ComplexF64}(undef, Nt)

        lattice = zeros(Int, (L,L)) 
        fname = "clock_"*"Nt=$mant"*"e$expo"*"L=$L"*"beta=$β"*".csv"
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
        end

        datafile = open(f1)
        data = DataFrame(
            E = E,
            m = m
        )
        CSV.write(f1, data)
        close(datafile)
        local elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
        println("$(round(now(), Dates.Second));\nβ = $β,$i/$(size(betarray,1)) elapsed time $(elapsed)")
    end

    elapsed = Dates.canonicalize(Dates.round((now()-start_all), Dates.Second))
    println("$(round(now(), Dates.Second))\nFINISHED L = $L\nTotal elapsed time = $(elapsed)")
end
