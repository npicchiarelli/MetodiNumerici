using Combinatorics, LinearAlgebra, Random


function init_lattice(L::Int ,q = 4)
    lattice = rand(0:q-1, (L,L))
    return lattice
end

function complex_vars(lattice::Array{Int}, q::Int = 4)
    return exp.(im.*2π.*lattice./q)
end

function getneighbors(idx::CartesianIndex, L::Int)
    lrud = [[0 1]; [0 -1]; [1 0]; [-1 0]] 
    i_v = hcat(getindex(idx,1), getindex(idx,2)) #matrice 1x2 con riga e colonna di idx
    nn_v = i_v .+ lrud # matrice 4x2 con riga e colonna del nearest neighbor di idx in ogni direzione
    replace!(nn_v, 0=>L, L+1=>1)#impongo le condizioni al bordo periodiche
    return CartesianIndex.(nn_v[:,1], nn_v[:,2])#ho un array con gli indici cartesiani dei nearest neig. di idx
end

function energy(lattice::Array{Int}, q::Int)
    e_persite = zeros(Float64, size(lattice))
    L = size(lattice, 1)
    for i in CartesianIndices(lattice)
        nn = getneighbors(i,L)
        e_persite[i] = -sum(cos, (2pi/q).*(lattice[i] .- lattice[nn]))
    end
    e_tot = sum(e_persite)
    return e_tot/L^2
end

function magnetization(lattice::Array{Int}, q::Int, L::Int)
    comp_lat = complex_vars(lattice, q)
    m = sum(comp_lat)/L^2
    return m
end

function energy_local(S::Int, lattice::Array{Int}, rx::Int, ry::Int, q::Int)
    L = size(lattice,1)
    nn = getneighbors(CartesianIndex(rx,ry), L)
    energy = -sum(cos, (2pi/q).*(S .- lattice[nn]))
    return energy
end

function init_prob_dict(q::Int, β::Float64)
    a = 0:q-1
    comb_arr = collect(with_replacement_combinations(a,4))#array con tutte le possibili combinazioni con ripetizione di q elementi e lunghezza 4
    e_v = Vector{Vector{Float64}}(undef, 35)#vettore di vettori di lunghezza bohx35 (35 è il numero di combinazioni)
    prob_v = Vector{Vector{Float64}}(undef, 35)
    norm_prob = Vector{Vector{Float64}}(undef, 35)

    for (idx, i) in enumerate(comb_arr)#i è l'elemento di comb_arr idx l'indice di comb_arr
        e_v[idx] = round.([-sum(cos, (2pi/q).*(j .- i)) for j in a])#= all'elemento e[idx] associo un array con il valore dell'energia per ogni 
        valore di j (valore centrale) è un array con j valori 
        =#
        prob_v[idx] = exp.(-β* e_v[idx])./sum(exp.(-β.*e_v[idx]))#calcolo probabilità associata ad ogni valore di j
        # @show i, e_v[idx]
    end
    
    for i in 1:length(prob_v)
        cum_prob[i] = ([sum(prob_v[i][1:j]) for j in 1:4])#calcola la somma delle probabilità da 0 a i e ho un array con 4 elementi

    return Dict(zip(comb_arr, cum_prob))#associo alla probabilità cumulativa la combinazione di nnb corrispondente
end

function metropolis!(lattice::Array{Int}, rx::Int, ry::Int, Δ::Int, q::Int, β::Float64)
    ΔS = Int(floor(rand(Float64)*Δ))
    S_test = (lattice[rx,ry] + ΔS)%q
    
    ΔE = energy_local(S_test, lattice, rx, ry,q) - energy_local(lattice[rx,ry], lattice, rx, ry,q) 

    if ΔE < 0
        lattice[rx,ry] = S_test
        return 1
    elseif rand(Float64) < exp(-β*ΔE)
        lattice[rx,ry] = S_test
        return 1
    end

    return 0
end

function heathbath!(lattice::Array{Int}, idx::CartesianIndex, pdict::Dict, L)
    acc = 0
    nn_idx = getneighbors(idx, L)
    nn = lattice[nn_idx]
    prob = pdict[sort(nn)]

    rnumber = rand()
    rnumber
    for (i,p) in enumerate(prob)
        if rnumber<p # se r è sotto al primo elemento metto 0, se è sotto il secondo metto 1 e così via
            lattice[idx] = i-1
            acc+=1
            break
        end
    end 
    return acc
end

function JackKnife(e::Array, m::Array, blocksize::Int, L::Int)
    length = size(e,1) ÷ blocksize 
    V = L*L

    e_cut = e[1:length*blocksize]
    m_cut = m[1:length*blocksize]

    e_b = reshape(e_cut, (blocksize, length))
    m_b = reshape(m_cut, (blocksize, length))

    e_blocked = mean(e_b, dims = 1)
    e2_blocked = mean(x->x^2,e_b, dims = 1)
    m_abs_blocked = mean(abs, m_b, dims = 1)
    m_abs2_blocked = mean(abs2, m_b, dims = 1)
    m_abs4_blocked = mean(x->abs(x)^4, m_b, dims = 1)

    jk(a) = let s = sum(a); [s-v for v in a] end

    e_j = vec(jk(e_blocked)/(length-1))
    e2_j = vec(jk(e2_blocked)/(length-1))
    m_abs_j = vec(jk(m_abs_blocked)/(length-1))
    m_abs2_j = vec(jk(m_abs2_blocked)/(length-1))
    m_abs4_j = vec(jk(m_abs4_blocked)/(length-1))

    
    spec_heat = V*(e2_j.-e_j.^2 )
    susc = V*(m_abs2_j.-m_abs_j.^2)
    bind = m_abs4_j./(m_abs2_j.^2)
    
    vars = [e_j m_abs_j V*m_abs2_j spec_heat susc bind]
    media = mean(vars, dims = 1)
    return media, stdm(vars, media, dims = 1, corrected = false).*sqrt(length-1)
end

function Blocking(e::Array, blocksize::Int)
    length = size(e,1) ÷ blocksize 
    e_cut = e[1:length*blocksize]
    e_b = reshape(e_cut, (blocksize, length))

    return mean(e_b, dims = 1), length
end