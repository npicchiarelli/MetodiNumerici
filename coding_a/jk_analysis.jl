using CSV, DataFrames, Dates, Plots, Statistics
include("clock.jl")

L = 60
β = 0.85
path = "..\\simulations_a\\" 

v = []
blocksize_v = collect(Int64, 100:100:10000)
# blocksize_v = 100

f1 = path*"clock_Nt=1e6"*"L=$L"*"beta=$β.csv"

for blocksize in blocksize_v
    df= CSV.read(f1, DataFrame, types = [Float64, ComplexF64])

    energ = df[!,:E]
    magn = df[!,:m]

    means, stds = JackKnife(energ,magn,blocksize)

    # eb, l = Blocking(energ, blocksize)
    push!(v,stds[1])
    println(blocksize)

end
p = plot()
xticks!(blocksize_v[begin:5:end])
plot!(blocksize_v, v)
