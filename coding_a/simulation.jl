using CSV, DataFrames
include("clock.jl")

L= 3
Q = 4
Nt = 1000
β = 1.
Nupdates = Nt*L*L

pdict = init_prob_dict(Q,β)
E = Vector{Float64}(undef, Nt)
lattice = zeros(Int, (L,L))

path = "..\\simulations_a\\" 
f1 = path*"clock.csv"
touch(f1)
acc = 0

for nt in 1:Nt
    for idx in CartesianIndices(lattice)
        global acc += heathbath!(lattice, idx, pdict, L) 
    end
E[nt] = energy(lattice, Q, L)
datafile = open(f1)
data = DataFrame(
    Energy = E
)
CSV.write(f1, data)
end