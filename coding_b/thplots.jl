using CSV, DataFrames, Dates, Plots, Printf, Statistics, LaTeXStrings
include("harm_osc.jl")

default(fontfamily = "Computer Modern",
background_color = :transparent,
foreground_color = :black,
margin=5Plots.mm
)

meanH(α::Float64) = 0.5*coth(α/2)
 
path = "..\\simulations_b\\"
sample = 100000000
fname = @sprintf "data%.1e.txt" sample

df= CSV.read(joinpath([path, fname]), DataFrame)

xp = plot(legend =:bottomright)
x2p = plot(legend =:right)
Kp = plot(legend = false)
Hp = plot(legend =:right)

cst = meanH(5.)

scatter!(xp, df[!,:eta_list].^2, df[!,:x], yerr = df[!,:xv], markershape=:plus, label = "Sampled")
plot!(xp, [0, 1.6], [0, 0], lw=1, lc=:black, ls =:dash, label = "Analytical")
title!(xp, "Average Position, βħω = 5.0")
ylabel!(xp, L"$\langle x \rangle$")
xlabel!(xp, L"$\eta^2$")

scatter!(x2p, df[!,:eta_list].^2, df[!,:x2], yerr = df[!,:x2v], markershape=:plus, label = "Sampled")
plot!(x2p, [0, 1.6], [cst, cst], lw=1, lc=:black, ls =:dash, label = L"\frac{1}{2}coth(\beta \hbar \omega /2)",)
title!(x2p, "Average Square Position, βħω = 5.0")
ylabel!(x2p, L"$\langle x^2 \rangle$")
xlabel!(x2p, L"$\eta^2$")

scatter!(Kp, df[!,:eta_list].^2, df[!,:K], yerr = df[!,:Kv], markershape=:plus)
title!(Kp, "Average (wrong) Kinetic Energy, βħω = 5.0")
ylabel!(Kp, L"$\langle K \rangle$")
xlabel!(Kp, L"$\eta^2$")

scatter!(Hp, df[!,:eta_list].^2, df[!,:H], yerr = df[!,:Hv], markershape=:plus, label = "Sampled")
plot!(Hp, [0, 1.6], [cst, cst], lw=1, lc=:black, ls =:dash, label = L"\frac{1}{2}coth(\beta \hbar \omega /2)")
title!(Hp, "Average Energy, βħω = 5.0")
ylabel!(Hp, L"$\langle H \rangle$")
xlabel!(Hp, L"$\eta^2$")

display(xp)
display(x2p)
display(Kp)
display(Hp)

# savefig(xp, "..\\imgs_b/xplot.png")
# savefig(x2p, "..\\imgs_b/x2plot.png")
# savefig(Hp, "..\\imgs_b/Hplot.png")
# savefig(Kp, "..\\imgs_b/Kplot.png")



