using CSV, DataFrames, Dates, Plots, Statistics

function JackKnife(e::Array, m::Array, blocksize::Int)
    length = size(e,1) ÷ blocksize 

    m_j = Vector{Float64}(undef, length) 
    e_j = Vector{Float64}(undef, length) 
    m_abs_j = Vector{Float64}(undef, length)
    m_abs2_j = Vector{Float64}(undef, length)
    m_abs4_j = Vector{Float64}(undef, length)
    e2_j = Vector{Float64}(undef, length)

    e_cut = e[1:length*blocksize]
    m_cut = m[1:length*blocksize]

    e_b = reshape(e_cut, (blocksize, length))
    m_b = reshape(m_cut, (blocksize, length))

    e_j = mean(e_b, dims = 1)
    e2_j = mean(x->x^2,e_b, dims = 1)
    m_abs_j = mean(abs, m_b, dims = 1)
    m_abs2_j = mean(abs2, m_b, dims = 1)
    m_abs4_j = mean(x->abs(x)^4, m_b, dims = 1)

    spec_heat = e2_j.-e_j.^2 
    susc = m_abs2_j.-m_abs_j.^2
    bind = m_abs4_j./(m_abs2_j.^2)
    
    vars = [e_j; m_abs_j; m_abs2_j; spec_heat; susc; bind]
    media = mean(vars, dims = 2)
    return media, stdm(vars, media, dims = 2).*sqrt(length)
end

L = 80
betarray = LinRange(0.85,0.91,21)
path = "..\\simulations_a\\" 
f_w =  path*"data"*"L=$L"*".csv"
touch(f_w)

tosave = DataFrame([[],[],[],[],[],[],[],[],[],[],[],[],[],[]], ["beta", "m", "e","e_v","m_abs","m_abs_v","m_abs2","m_abs2_v","spec_heat","spec_heat_v","susc","susc_v","bind","bind_v"])

for (i,β) in enumerate(betarray)
    # start = now()
    local f1 = path*"clock_Nt=1e5"*"L=$L"*"beta=$β.csv"
    df= CSV.read(f1, DataFrame, types = [Float64, ComplexF64])
    # @show df
    # elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
    # println("Opening took: $(elapsed)")

    energ = df[!,:E]
    magn = df[!,:m]

    start = now()

    @show means, stds = JackKnife(energ,magn,1000)

    inter = ([means stds]'[:])'
    data = [β abs(mean(magn)) inter]

    push!(tosave, data)
    println("β = $β,$i/$(size(betarray,1))")
end


CSV.write(f_w, tosave)