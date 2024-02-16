using CSV, DataFrames, Dates
include("clock.jl")

L= 4
Q = 4
expo = 7
Nt = 1*Int(10^expo)
#β = 0.42
Δ = 1

metropolis = false

@show Nupdates = Nt*L*L
for β in []
    pdict = init_prob_dict(Q,β)
    E = Vector{Float64}(undef, Nt)
    m = Vector{ComplexF64}(undef, Nt)

    lattice = zeros(Int, (L,L))

    path = "..\\simulations_a\\" 
    datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")

    f1 = path*"clock_"*"Nt=1e$expo"*"L=$L"*"beta=$β"*".csv"
    # f1 = path*"clock_"*"L=$L"*"Nt=1e$expo"*"_$datestamp"*".csv"

    touch(f1)
    acc = 0
    start = now()
    for nt in 1:Nt
        # start_step = now()
        for idx in CartesianIndices(lattice)
            if metropolis
                rx = idx[1]
                ry = idx[2]
                global acc += metropolis!(lattice,rx,ry,Δ,L,β) 
            else
                global acc += heathbath!(lattice, idx, pdict, L)
            end
        end
        if nt%100000 == 0
            elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
            acc_p = acc/(nt*L*L)*100
            # per_step = Dates.canonicalize(now()-start_step)
            println("Step $nt, acceptance = $acc_p %, total elapsed time $(elapsed)")#, time per step $per_step")
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
end
# pdict = init_prob_dict(Q,β)
# E = Vector{Float64}(undef, Nt)
# m = Vector{ComplexF64}(undef, Nt)

# lattice = zeros(Int, (L,L))

# path = "..\\simulations_a\\" 
# datestamp=Dates.format(now(),"YYYYmmdd-HHMMSS")

# f1 = path*"clock_"*"Nt=1e$expo"*"L=$L"*"beta=$β"*".csv"
# # f1 = path*"clock_"*"L=$L"*"Nt=1e$expo"*"_$datestamp"*".csv"

# touch(f1)
# acc = 0
# start = now()
# for nt in 1:Nt
#     # start_step = now()
#     for idx in CartesianIndices(lattice)
#         if metropolis
#             rx = idx[1]
#             ry = idx[2]
#             global acc += metropolis!(lattice,rx,ry,Δ,L,β) 
#         else
#             global acc += heathbath!(lattice, idx, pdict, L)
#         end
#     end
#     if nt%100000 == 0
#         elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
#         acc_p = acc/(nt*L*L)*100
#         # per_step = Dates.canonicalize(now()-start_step)
#         println("Step $nt, acceptance = $acc_p %, total elapsed time $(elapsed)")#, time per step $per_step")
#     end

# E[nt] = energy(lattice, Q, L)
# m[nt] = magnetization(lattice, Q, L)
# end

# datafile = open(f1)
# data = DataFrame(
#     E = E,
#     m=m
# )
# CSV.write(f1, data)

