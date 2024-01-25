using Combinatorics
q = 4
a = 0:q-1
comb_arr = collect(with_replacement_combinations(a,4))

for i in comb_arr
    e_v = [-sum(cos, (2pi/q).*(j .- i)) for j in a]
    @show i, e_v
end