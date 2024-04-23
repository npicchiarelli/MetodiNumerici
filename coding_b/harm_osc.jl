using LinearAlgebra

function calc_x(lattice::Array{Float64}, Nt::Int)
    return sum(lattice)/Nt    
end

function calc_x2(lattice::Array{Float64}, Nt::Int)
    return sum(abs2,lattice)/Nt
end


function calc_Knaive(lattice::Array{Float64}, Nt::Int, eta::Float64)
    diffs = lattice .- circshift(lattice, 1)
    return sum(diffs.*diffs)/(2*Nt*eta*eta)
end

function metropolis!(lattice::Array{Float64}, r::Int, Δ::Int, η::Float64)
    ΔS = Δ*(2*rand(Float64)-1)
    trial = (lattice[r]+ΔS)
    Eold = lattice[r]*lattice[r]*(η/2.0+1.0/η)-lattice[r]*(circshift(lattice, 1)[r]+circshift(lattice, -1)[r])/η
    Enew = trial*trial*(η/2.0+1.0/η)-trial*(circshift(lattice, 1)[r]+circshift(lattice, -1)[r])/η

    if Enew < Eold
        lattice[r] = trial
        return 1
    elseif rand(Float64) < exp(-(Enew-Eold))
        lattice[r] = trial
        return 1
    end

    return 0
end

function heathbath!(lattice::Array{Float64}, r::Int, eta::Float64)
    std = 1.0/sqrt(eta + 2.0/eta)
    avg = (circshift(lattice, 1)[r]+circshift(lattice, -1)[r])/(eta*(eta + 2.0/eta))
    
    lattice[r] = avg+std*randn()

    return 1
end

function overrelax!(lattice::Array{Float64}, r::Int, eta::Float64)
    avg = (circshift(lattice, 1)[r]+circshift(lattice, -1)[r])/(eta*(eta + 2.0/eta))
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
