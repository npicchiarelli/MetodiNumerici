using Combinatorics, LinearAlgebra, Random

Random.seed!(666)
function init_lattice(L::Int ,q = 4)
    lattice = rand(0:q-1, (L,L))
    return lattice
end
function getneighbors(idx::CartesianIndex, L::Int)
    lrud = [[0 1]; [0 -1]; [1 0]; [-1 0]] 
    i_v = hcat(getindex(idx,1), getindex(idx,2))
    nn_v = i_v .+ lrud
    replace!(nn_v, 0=>L, L+1=>1)
    return CartesianIndex.(nn_v[:,1], nn_v[:,2])
end

function energy(lattice::Array{Int}, q::Int)
    e_persite = zeros(Float64, size(lattice))
    L = size(lattice, 1)
    for i in CartesianIndices(lattice)
        nn = getneighbors(i,L)
        e_persite[i] = -sum(cos, (2pi/q).*(lattice[i] .- lattice[nn]))
    end
    e_tot = sum(e_persite)
    return e_tot
end

function energy_local(S::Int, lattice::Array{Int}, rx::Int, ry::Int, q::Int)
    L = size(lattice,1)
    nn = getneighbors(CartesianIndex(rx,ry), L)
    energy = -sum(cos, (2pi/q).*(S .- lattice[nn]))
    return energy
end

function init_prob_dict(q::Int, β::Float64)
    comb_arr = collect(with_replacement_combinations(a,4))
    e_v = Vector{Vector{Float64}}(undef, 35)
    prob_v = Vector{Vector{Float64}}(undef, 35)
    for (idx, i) in enumerate(comb_arr)
        e_v[idx] = round.([-sum(cos, (2pi/q).*(j .- i)) for j in a])
        prob_v[idx] = exp.(-e_v[idx])./sum(exp.(-e_v[idx]))
        # @show i, e_v[idx]
    end

    return Dict(zip(comb_arr, prob_v))
end

function metropolis!(lattice::Array{Int}, rx::Int, ry::Int, Δ::Int, q::Int)
    acc = 0
    @show ΔS = Int(floor(rand(Float64)*Δ))
    S_test = mod((lattice[rx,ry] + ΔS), q)
    
    ΔE = energy_local(lattice[rx,ry], lattice, rx, ry,q) - energy_local(S_test, lattice, rx, ry,q) 

    if ΔE < 0
        lattice[rx,ry] = S_test
        acc+=1
    elseif rand(Float64) < exp(ΔE)
        lattice[rx,ry] = S_test
        acc+=1
    end

    return acc
end

function heathbath!(lattice::Array{Int}, idx::CartesianIndex, q::Int, pdict::Dict, L)
    nn_idx = getneighbors(idx, L)
    nn = lattice[nn_idx]
    prob = pdict[sort(nn)]

    rnumber = rand(Float64)
    for (i,p) in enumerate(prob)
        if rnumber<p
            lattice[idx] = i-1
            break
    end 
    return 1
end

function findedict(lattice::Array{Int}, rx::Int, ry::Int, edict::Dict)
    L = size(lattice, 1)
    lrud = [[0 1]; [0 -1]; [1 0]; [-1 0]] 
    i_v = [rx ry]
    nn_v = i_v .+ lrud
    replace!(nn_v, 0=>L, L+1=>1)
    nn = CartesianIndex.(nn_v[:,1], nn_v[:,2])

    neighbors = sort(lattice[nn])
    return edict[neighbors]
end

Q = 4
L = 3
Δ = 4

@show lattice = init_lattice(L)
energy_dict = init_energy_dict(Q, 1.)


# findedict(lattice, 2, 2, energy_dict)



#=
@show energy(lattice, 4)
rx = rand(1:3)
ry = rand(1:3)
println(rx,ry)
@show metropolis!(lattice, rx, ry, Δ, 4)
@show lattice
@show energy(lattice, 4)=#