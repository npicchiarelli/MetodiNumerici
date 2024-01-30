using Combinatorics
q = 4
a = 0:q-1
comb_arr = collect(with_replacement_combinations(a,4))
e_v = Vector{Vector{Float64}}(undef, 35)

for (idx, i) in enumerate(comb_arr)
    e_v[idx] = round.([-sum(cos, (2pi/q).*(j .- i)) for j in a])
    @show i, e_v[idx]
end

# energy_dict = Dict(zip(comb_arr, e_v))