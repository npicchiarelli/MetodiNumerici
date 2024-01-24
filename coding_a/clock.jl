using LinearAlgebra, Random

# Random.seed!(42)
function init_lattice(L::Int ,q = 4)
    lattice = rand(0:q-1, (L,L))
    return lattice
end

function energy(lattice::Array{Int}, q::Int)
    e_persite = zeros(Float64, size(lattice))
    lrud = [[0 1]; [0 -1]; [1 0]; [-1 0]] 
    L = size(lattice, 1)
    for i in CartesianIndices(lattice)
        i_v = hcat(getindex(i,1) ,getindex(i,2))
        nn_v = i_v .+ lrud
        replace!(nn_v, 0=>L, L+1=>1)
        nn = CartesianIndex.(nn_v[:,1], nn_v[:,2])
        e_persite[i] = sum(lattice[i] .- lattice[nn])
        # e_persite[i] = sum(cos, (2pi/q).*(lattice[i] .- lattice[nn]))
    end
    return e_persite
end

@show lattice = init_lattice(2)
energy(lattice, 4)

