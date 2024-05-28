using CSV, DataFrames, Dates, Plots, Printf, Statistics, LaTeXStrings
include("harm_osc.jl")

default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
margin=5Plots.mm
)

meanH(α::Float64) = 0.5*coth(α/2) 
 
path = "..\\simulations_b\\"
sample = 100000000
fname = @sprintf "data%.1e.txt" sample

df= CSV.read(joinpath([path, fname]), DataFrame)

xp = plot(legend = false)
x2p = plot(legend = false)
Kp = plot(legend = false)
Hp = plot(legend = false)

cst = meanH(5.)

scatter!(xp, df[!,:eta_list].^2, df[!,:x], yerr = df[!,:xv], markershape=:plus)
plot!(xp, [0, 1.6], [0, 0], lw=1, lc=:black, ls =:dash, label = false)
title!(xp, "Average Position")
ylabel!(xp, L"$\langle x \rangle$")
xlabel!(xp, L"$\eta^2$")

scatter!(x2p, df[!,:eta_list].^2, df[!,:x2], yerr = df[!,:x2v], markershape=:plus)
plot!(x2p, [0, 1.6], [cst, cst], lw=1, lc=:black, ls =:dash, label = false)
title!(x2p, "Average Square Position")
ylabel!(x2p, L"$\langle x^2 \rangle$")
xlabel!(x2p, L"$\eta^2$")

scatter!(Kp, df[!,:eta_list].^2, df[!,:K], yerr = df[!,:Kv], markershape=:plus)
title!(Kp, "Average (wrong) Kinetic Energy")
ylabel!(Kp, L"$\langle K \rangle$")
xlabel!(Kp, L"$\eta^2$")

scatter!(Hp, df[!,:eta_list].^2, df[!,:H], yerr = df[!,:Hv], markershape=:plus)
plot!(Hp, [0, 1.6], [cst, cst], lw=1, lc=:black, ls =:dash, label = false)
title!(Hp, "Average Energy")
ylabel!(Hp, L"$\langle H \rangle$")
xlabel!(Hp, L"$\eta^2$")

display(xp)
display(x2p)
display(Kp)
display(Hp)

