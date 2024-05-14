using DelimitedFiles, Plots, Printf, Statistics
include("harm_osc.jl")

path = "..\\simulations_b\\"
sample, Nt, β = 100000000, 60, 10.
fname = @sprintf "data_sample%.1eNt%ibeta%.1f.csv" sample Nt β

df= CSV.read(joinpath([path, fname]), DataFrame)

p = plot()

scatter!(p, df[!,:sp],df[!,:x], yerr = df[!,:xv], markershape=:plus)
scatter!(p, df[!,:sp],df[!,:x2], yerr = df[!,:x2v], markershape=:plus)
scatter!(p, df[!,:sp],df[!,:x3], yerr = df[!,:x3v], markershape=:plus)
scatter!(p, df[!,:sp],df[!,:x4], yerr = df[!,:x4v], markershape=:plus)

display(p)