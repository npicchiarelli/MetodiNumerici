using CSV, DataFrames,Plots, Statistics

function blocking(x::Array, blocksize::Int)
    length = size(x,1) ÷ blocksize 
    x_b = Vector{Float64}(undef, length)
    x_cut = x[1:length*blocksize]
    x_b = mean(reshape(x_cut, (blocksize, length)), dims = 1)
    
    return x_b
end

function JackKnife(e::Array, m::Array, blocksize::Int)
    length = size(e,1) ÷ blocksize 
    m_abs = Vector{Float64}(undef, length)
    m_abs2 = Vector{Float64}(undef, length)
    m_abs4 = Vector{Float64}(undef, length)
    e2 = Vector{Float64}(undef, length)

    e_cut = e[1:length*blocksize]
    m_cut = m[1:length*blocksize]

    e_b = reshape(e_cut, (blocksize, length))
    m_b = reshape(m_cut, (blocksize, length))

    e_j = [mean([e[1:i*blocksize]; e[(i+1)*blocksize:end]]) for i = 0:length-1]
    e2_j = [mean(abs2, [e[1:i*blocksize]; e[(i+1)*blocksize:end]]) for i = 0:length-1]
    m_abs_j = [mean(abs, [m[1:i*blocksize]; m[(i+1)*blocksize:end]]) for i = 0:length-1]
    m_abs2_j = [mean(abs2, [m[1:i*blocksize]; m[(i+1)*blocksize:end]]) for i = 0:length-1]
    m_abs4_j = [mean((abs.([m[1:i*blocksize]; m[(i+1)*blocksize:end]])).^4) for i = 0:length-1]

    spec_heat = e2_j.-e_j.^2 
    susc = m_abs2_j.-m_abs_j.^2
    bind = m_abs4_j./(m_abs2_j.^2)
    
    vars = [e_j m_abs_j m_abs2_j spec_heat susc bind]
    media = mean(vars, dims = 1)
    return media, stdm(vars, media, dims = 1)*sqrt(length)
end




path = "..\\simulations_a\\" 
f1 = path*"clock_L=4Nt=1e5_20240201-111634.csv"
L = 4
df= CSV.read(f1, DataFrame, types = [Float64, ComplexF64])
# @show df

e = df[!,:E]
m = df[!,:m]

@show mean(m)

JackKnife(e,m,100)