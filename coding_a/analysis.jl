using CSV, DataFrames, Dates, Plots, Statistics
include("clock.jl")

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

    @show means, stds = JackKnife(energ,magn,100)

    inter = ([means stds]'[:])'
    data = [β abs(mean(magn)) inter]

    push!(tosave, data)
    println("β = $β,$i/$(size(betarray,1))")
end


CSV.write(f_w, tosave)