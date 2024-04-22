function calc_x(lattice::Array{Float64}, Nt::Int)
    return sum(lattice)/Nt    
end

function calc_x2(lattice::Array{Float64}, Nt::Int)
    return sum(abs2,lattice)/Nt
end

function init_neighbors(Nt::Int)
    idx = collect(range(1,Nt,Nt))
    nnl = idx.-1
    nnr = idx.+1
    replace!(nnl, 0=>Nt)
    replace!(nnr, Nt+1=>0)
    return Int.(nnl), Int.(nnr)
end

function calc_Knaive(lattice::Array{Float64}, Nt::Int, eta::Float64)
    diffs = lattice .- circshift(lattice, 1)
    return sum(diffs.*diffs)/(2*Nt*eta*eta)
end

function metropolis!(lattice::Array{Int}, r::IndexLinear, Δ::Int, β::Float64, η::Float64, nnl::Int, nnr::Int)
    ΔS = Δ*(2*rand(Float64)-1)
    trial = (lattice[r]+ΔS)
    Eold = lattice[r]*lattice[r]*(η/2.+1./η)-lattice[r]*(lattice[nnl[r]]+lattice[nnr[r]])/η
    Enew = trial*trial*(η/2.+1./η)-trial*(lattice[nnl[r]]+lattice[nnr[r]])/η

    if Enew < Eold
        lattice[r] = trial
        return 1
    elseif rand(Float64) < exp(-(Enew-Eold))
        lattice[r] = trial
        return 1
    end

    return 0
end

function heathbath!(lattice::Array{Int}, r::IndexLinear, eta::Float64, nnl::Int, nnr::Int)
    std = 1./sqrt(eta + 2./eta)
    avg = (lattice[nnl[r]]+lattice[nnr[r]])/(eta*(eta + 2./eta))
    
    lattice[r] = 
end
