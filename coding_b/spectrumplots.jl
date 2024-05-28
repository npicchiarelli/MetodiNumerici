using CSV, DataFrames, DelimitedFiles, Plots, Printf, Statistics
include("harm_osc.jl")

default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
background_color_legend = nothing,
margin=5Plots.mm
)

path = "..\\simulations_b\\"
sample, Nt, β = 100000000, 60, 10.
fname = @sprintf "data_sample%.1eNt%ibeta%.1f.csv" sample Nt β

df= CSV.read(joinpath([path, fname]), DataFrame)

p = plot(legend=:best)

scatter!(p, df[!,:sp],df[!,:x], yerr = df[!,:xv], markershape=:plus, label = L"x-x")
scatter!(p, df[!,:sp],df[!,:x2], yerr = df[!,:x2v], markershape=:plus, label = L"x^2-x^2")
scatter!(p, df[!,:sp],df[!,:x3], yerr = df[!,:x3v], markershape=:plus, label = L"x^3-x^3")
scatter!(p, df[!,:sp],df[!,:x4], yerr = df[!,:x4v], markershape=:plus, label = L"A-A")
xlabel!("Spacing")
ylabel!("Correlators")
ylims!(0.9, 4)
annotate!(p, 2.8, 3.33, text(L"Where $A = x^3-\frac{3}{2}x$", :black,:left, 8))
title!("Plot of Correlators vs Spacing")

display(p)