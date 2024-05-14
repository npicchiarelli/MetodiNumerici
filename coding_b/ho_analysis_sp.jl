using DelimitedFiles, Plots, Printf, Statistics
include("harm_osc.jl")

path = "..\\simulations_b\\old\\"
blocksize = 1000
therm = 10000
sample = 100000000
startp = @sprintf "ho_sp_sample=%.1e" sample
paths = filter(startswith(startp), readdir(path))

for fname in [paths[2]]
    bend = findfirst("beta=", fname)[end]
    Nstr = findfirst("Nt=", fname)[1]
    Nend = findfirst("Nt=", fname)[end]
    Nt = parse(Int, fname[Nend+1:end-4])
    β = parse(Float64, fname[bend+1:Nstr-1])
    η = β/Nt

    local w = open(joinpath([path, fname]), "r") do io
        readdlm(io, header = true)
    end
    
    spacing=Vector{Int}()
    for el in w[2][2:end]
        append!(spacing, parse(Int, el[4:end]))
    end

    unique!(spacing)
    cormat = w[1][therm+1:end,:]
    datajack = zeros(size(cormat,1)÷blocksize, size(cormat,2))

    for i in 1:size(cormat,2)
        datajack[:,i] = JackKnife(cormat[:,i], blocksize)
    end

    corrected = datajack
    corrected[:,3:4:end].-=(corrected[:,1].^2)
    corrected = corrected[:,2:end]
    gaps = log.(corrected./circshift(corrected, (0,-4)))
    errs = std(gaps, dims = 1, corrected = false).*sqrt(size(datajack,1)-1)
    gaps = mean(gaps, dims = 1)
    gaps = gaps[:,1:end-4]./η
    errs = errs[:,1:end-4]./η

    # p = plot()
    # scatter!(p,spacing[1:end-1], gaps[1:4:end], yerr = errs[1:4:end], markershape=:plus)
    # scatter!(p,spacing[1:end-1], gaps[2:4:end], yerr = errs[2:4:end], markershape=:plus)
    # scatter!(p,spacing[1:end-1], gaps[3:4:end], yerr = errs[3:4:end], markershape=:plus)
    # scatter!(p,spacing[1:end-1], gaps[4:4:end], yerr = errs[4:4:end], markershape=:plus)
    # display(p)
end