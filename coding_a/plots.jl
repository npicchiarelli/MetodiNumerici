using CSV, DataFrames, Dates, Plots, Statistics

p = plot()
for size in [20 40 50 60 70]
    f = "..\\simulations_a\\data\\dataL=$size.csv"

    df= CSV.read(f, DataFrame)
    println(names(df))

    x = df[!,:beta]
    y = df[!,:m_abs2]
    sy = df[!,:m_abs2_v]
    plot!(p,x,y, yerr = sy, label = "$size")
end

p