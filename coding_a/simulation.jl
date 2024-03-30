using CSV, DataFrames, Dates
include("clock.jl")

L= 80
Q = 4
expo = 5
Nt = 1*Int(10^expo)
#β = 0.42
Δ = 1
path = "..\\simulations_a\\"
metropolis = false

@show Nupdates = Nt*L*L

betarray = round.(LinRange(0.85,0.91,41), digits = 4)
start_all = now()
for (i, β) in enumerate(betarray[34:end])
    pdict = init_prob_dict(Q, β)
    E = Vector{Float64}(undef, Nt)
    m = Vector{ComplexF64}(undef, Nt)

    lattice = zeros(Int, (L,L)) 
    # datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")

    f1 = path*"clock_"*"Nt=1e$expo"*"L=$L"*"beta=$β"*".csv"

    touch(f1)
    acc = 0
    start = now()
    for nt in 1:Nt
        # start_step = now()
        for idx in CartesianIndices(lattice)
            if metropolis
                rx = idx[1]
                ry = idx[2]
                acc += metropolis!(lattice,rx,ry,Δ,L,β) 
            else
                acc += heathbath!(lattice, idx, pdict, L)
            end
        end
        # if nt%100000 == 0
        #     elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
        #     println("Step $nt, acceptance = $acc_p %, total elapsed time $(elapsed)")#, time per step $per_step")
        # end
    E[nt] = energy(lattice, Q, L)
    m[nt] = magnetization(lattice, Q, L)
    end

    datafile = open(f1)
    data = DataFrame(
        E = E,
        m = m
    )
    CSV.write(f1, data)
    local elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
    println("β = $β,$i/$(size(betarray,1)) elapsed time $(elapsed)")
end

elapsed = Dates.canonicalize(Dates.round((now()-start_all), Dates.Second))
println("FINISHED\nTotal elapsed time = $(elapsed)")
