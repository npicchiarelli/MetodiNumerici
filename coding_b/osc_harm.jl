using LinearAlgebra

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

function metropolis!(lattice::Array{Float64}, r::IndexLinear, Δ::Int, β::Float64, η::Float64, nnl::Int, nnr::Int)
    ΔS = Δ*(2*rand(Float64)-1)
    trial = (lattice[r]+ΔS)
    Eold = lattice[r]*lattice[r]*(η/2.0+1.0/η)-lattice[r]*(lattice[nnl[r]]+lattice[nnr[r]])/η
    Enew = trial*trial*(η/2.0+1.0/η)-trial*(lattice[nnl[r]]+lattice[nnr[r]])/η

    if Enew < Eold
        lattice[r] = trial
        return 1
    elseif rand(Float64) < exp(-(Enew-Eold))
        lattice[r] = trial
        return 1
    end

    return 0
end

function heathbath!(lattice::Array{Float64}, r::IndexLinear, eta::Float64, nnl::Int, nnr::Int)
    std = 1.0/sqrt(eta + 2.0/eta)
    avg = (lattice[nnl[r]]+lattice[nnr[r]])/(eta*(eta + 2.0/eta))
    
    lattice[r] = avg+std*randn()

    return 1
end

function overrelax(lattice::Array{Float64}, r::IndexLinear, eta::Float64, nnl::Int, nnr::Int)
    avg = (lattice[nnl[r]]+lattice[nnr[r]])/(eta*(eta + 2.0/eta))
    new = 2.0*avg-lattice[r]
    lattice[r] = new

    return 1
end

function correlators(lattice::Array{Float64}, δt::Int)
    corr1 = dot(lattice, circshift(lattice, δt))
    corr2 = dot(lattice.^2,circshift(lattice, δt).^2)
    corr3 = dot(lattice.^3,circshift(lattice, δt).^3)
    corr3c = dot(map(x->x^3-1.5x, lattice), map(x->x^3-1.5x, circshift(lattice, δt)))
    return [corr1 corr2 corr3 corr3c]
end
