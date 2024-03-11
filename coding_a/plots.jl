using CSV, DataFrames, Dates, Plots, Statistics

p = plot()
for size in [80]
    f = "..\\simulations_a\\dataL=$size.csv"

    df= CSV.read(f, DataFrame)

    x = df[!,:beta]
    y = df[!,:spec_heat]
    sy = df[!,:spec_heat_v]
    plot!(p,x,y, yerr = sy, )
end

p