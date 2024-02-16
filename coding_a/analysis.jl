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

    e_b = mean(reshape(e_cut, (blocksize, length)), dims = 1)
    m_b = mean(reshape(m_cut, (blocksize, length)), dims = 1)

    for i = 0:length-1
        e_j[i+1] = mean([e_b[1:i]; e_b[(i+1):end]])
        e2_j[i+1] = mean(abs2, [e_b[1:i]; e_b[(i+1):end]])
        m_abs_j[i+1] = mean(abs, [m_b[1:i]; m_b[(i+1):end]])
        m_abs2_j[i+1] = mean(abs2, [m_b[1:i]; m_b[(i+1):end]])
        m_abs4_j[i+1] = mean((abs.([m_b[1:i]; m_b[(i+1):end]])).^4)
    end

    spec_heat = e2_j.-e_j.^2 
    susc = m_abs2_j.-m_abs_j.^2
    bind = m_abs4_j./(m_abs2_j.^2)
    
    vars = [e_j m_abs_j m_abs2_j spec_heat susc bind]
    media = mean(vars, dims = 1)
    return media, stdm(vars, media, dims = 1)*sqrt(length)
end

L = 4
betarray = LinRange(0.4,0.47,10)
path = "..\\simulations_a\\" 
f_w =  path*"data"*"L=$L"*".csv"
touch(f_w)

tosave = DataFrame([[],[],[],[],[],[],[],[],[],[],[],[],[],[]], ["beta", "m", "e","e_v","m_abs","m_abs_v","m_abs2","m_abs2_v","spec_heat","spec_heat_v","susc","susc_v","bind","bind_v"])

for β in betarray
    # start = now()
    f1 = path*"clock_Nt=1e4"*"L=$L"*"beta=$β.csv"
    df= CSV.read(f1, DataFrame, types = [Float64, ComplexF64])
    # @show df
    # elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
    # println("Opening took: $(elapsed)")

    energ = df[!,:E]
    magn = df[!,:m]

    # start = now()

    means, stds = JackKnife(energ,magn,1000)

    inter = ([means stds]'[:])'
    data = [β abs(mean(magn)) inter]

    # elapsed = Dates.canonicalize(Dates.round((now()-start), Dates.Second))
    # println("Computation took: $(elapsed)")
    # data = DataFrame(
    #     beta = β,
    #     m = mean(magn),
    #     e = means[1],
    #     e_v = stds[1],
    #     m_abs = means[2],
    #     m_abs_v = stds[2],
    #     m_abs2 = means[3],
    #     m_abs2_v = stds[3],
    #     spec_heat = means[4],
    #     spec_heat_v = stds[4],
    #     susc = means[5],
    #     susc_v = stds[5],
    #     bind = means[6],
    #     bind_v = stds[6],
    # )

    push!(tosave, data)
end

CSV.write(f_w, tosave)