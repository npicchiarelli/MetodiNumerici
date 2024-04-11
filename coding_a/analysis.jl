using CSV, DataFrames, Dates, Plots, Statistics
include("clock.jl")

size_list = [20, 40, 50, 60, 70, 80]
size_list = 20


for L in size_list
    if L == 60
        betarray = LinRange(0.85,0.91,41)
    else
        @show betarray =round.(LinRange(0.85,0.91,41), digits = 4)
    end
    path = "..\\simulations_a\\" 
    f_w =  path*"data"*"L=$L"*".csv"
    touch(f_w)

    tosave = DataFrame([[],[],[],[],[],[],[],[],[],[],[],[],[],[]], ["beta", "m", "e","e_v","m_abs","m_abs_v","m_abs2","m_abs2_v","spec_heat","spec_heat_v","susc","susc_v","bind","bind_v"])

    for (i,β) in enumerate(betarray[1:end])
        # start = now()
        local f1 = path*"clock_Nt=1e5"*"L=$L"*"beta=$β.csv"
        df= CSV.read(f1, DataFrame, types = [Float64, ComplexF64])

        energ = df[!,:E]
        magn = df[!,:m]

        therm = 10000
        start = now()

        
        means, stds = JackKnife(energ[therm:end],magn[therm:end],2000, L)
        inter = ([means' stds']'[:])'
        data = [β abs(mean(magn)) inter]

        push!(tosave, data)
        println("$L: β = $β,$i/$(size(betarray,1))")
    end


    CSV.write(f_w, tosave)
end


