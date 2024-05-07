using Statistics, DelimitedFiles, LinearAlgebra, Plots, Printf
include("harm_osc.jl")

therm = 1
blocksize = 100
path = joinpath(["..", "simulations_b"])
sample = 10000000
startp = @sprintf "ho_th_sample=%.1e" sample
paths = filter(startswith(startp), readdir(path))

size_list = []
eta_list = []
x = []
xv = []
x2 = []
x2v = []
K = []
Kv = []
H = []
Hv = []

for (i,fname) in enumerate(paths)
    Nt = parse(Int, fname[end-5:end-4])
    β = parse(Float64, fname[end-11:end-9])
    η = β/Nt

    w = open(joinpath([path, fname]), "r") do io
        readdlm(io, header = true)
    end

    x_j = JackKnife(w[1][therm:end,1], blocksize)
    x2_j = JackKnife(w[1][therm:end,2], blocksize)
    K_j = JackKnife(w[1][therm:end,3], blocksize)
    H_j = 0.5.*x2_j.-K_j .+ 1/(2*η)
    
    push!(size_list, Nt)
    push!(eta_list, η)
    push!(x, mean(x_j))
    push!(x2, mean(x2_j))
    push!(K, mean(K_j))
    push!(H, mean(H_j))

    push!(xv, std(x_j, corrected = false).*sqrt(length(x_j)-1))
    push!(x2v, std(x2_j, corrected = false).*sqrt(length(x_j)-1))
    push!(Kv, std(K_j, corrected = false).*sqrt(length(x_j)-1))
    push!(Hv, std(H_j, corrected = false).*sqrt(length(x_j)-1))
end

w = open(joinpath([path, "data.txt"]), "w") do io
    writedlm(io, ["size_list" "eta_list" "x" "xv" "x2" "x2v" "K" "Kv" "H" "Hv"], " ")
    writedlm(io, [size_list eta_list x xv x2 x2v K Kv H Hv], " ")
end


# scatter(eta_list.^2, H, yerr = Hv, size = (800, 500))
