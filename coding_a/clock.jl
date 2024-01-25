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
        # e_persite[i] = sum(lattice[i] .- lattice[nn])
        e_persite[i] = -sum(cos, (2pi/q).*(lattice[i] .- lattice[nn]))
    end
    e_tot = sum(e_persite)
    return e_tot
end

function energy_local(S::Int, lattice::Array{Int}, rx::Int, ry::Int, q::Int)
    i_v = [rx ry]
    lrud = [[0 1]; [0 -1]; [1 0]; [-1 0]]
    nn_v = i_v .+ lrud
    replace!(nn_v, 0=>L, L+1=>1)
    nn = CartesianIndex.(nn_v[:,1], nn_v[:,2])
    i = CartesianIndex(rx,ry)
    energy = -sum(cos, (2pi/q).*(S .- lattice[nn]))
    return energy
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

function heathbath!(lattice::Array{Int})
    
end

L = 3
Δ = 4
@show lattice = init_lattice(L)
@show energy(lattice, 4)
rx = rand(1:3)
ry = rand(1:3)
println(rx,ry)
@show metropolis!(lattice, rx, ry, Δ, 4)
@show lattice
@show energy(lattice, 4)