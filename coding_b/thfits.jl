using CSV, DataFrames, LsqFit, Plots
default(fontfamily = "Computer Modern",
background_color = :white,
foreground_color = :black,
margin=5Plots.mm
)

df = CSV.read("..\\simulations_b\\data1.0e+08.txt", DataFrame,)

esq = df[!,:eta_list].^2
fitotx2 = plot()
fitotH = plot()


# scatter!(xp,esq, df[!,:x], yerr = df[!,:xv])
# scatter!(x2p,esq, df[!,:x2], yerr = df[!,:x2v])
# scatter!(Kp,esq, df[!,:K], yerr = df[!,:Kv])
# scatter!(Hp,esq, df[!,:H], yerr = df[!,:Hv])

model(x,p) = p[1] .+ p[2].*x.^2

p2 = []
pH = []
p2v = []
pHv = []

etamax = [0.05, 0.1, 0.15, 0.2]

for maxeta in etamax
    msk = esq.<maxeta
    x2 = df[!,:x2][msk]
    x2v = df[!,:x2v][msk]
    H = df[!,:H][msk]
    Hv = df[!,:Hv][msk]
    eta = df[!,:eta_list][msk]

    fitx2 = curve_fit(model, eta, x2, x2v.^(-2), [0., 0.])
    println("x2")
    println(fitx2.param)
    covx2 = estimate_covar(fitx2)
    println("Errors on params: $(sqrt.(diag(covx2)))")


    fitH = curve_fit(model, eta, H, Hv.^(-2), [0., 0.])
    println("H")
    println(fitH.param)
    covH = estimate_covar(fitH)
    println("Errors on params: $(sqrt.(diag(covH)))")
    push!(p2, fitx2.param[1])
    push!(p2v, sqrt(covx2[1]))
    push!(pH, fitH.param[1])
    push!(pHv, sqrt(covH[1]))

end
meanH(α::Float64) = 0.5*coth(α/2) 
cst = meanH(5.)

scatter!(fitotx2, etamax, p2, yerr = p2v, markershape = :plus, label = "Sampled")
plot!(fitotx2, [0, 0.25], [cst, cst], lw=1, lc=:black, ls =:dash, label = L"\frac{1}{2}coth(\beta \hbar \omega /2)")
title!(fitotx2, L"$\langle x^2 \rangle$ extracted by fit, βħω = 5.0")
ylabel!(fitotx2, L"$\langle x^2 \rangle$")
xlabel!(fitotx2, L"$\eta^2_{max}$")

scatter!(fitotH, etamax, pH, yerr = pHv, markershape = :plus, label = "Sampled")
plot!(fitotH, [0, 0.25], [cst, cst], lw=1, lc=:black, ls =:dash, label = L"\frac{1}{2}coth(\beta \hbar \omega /2)")
title!(fitotH, L"$\langle H \rangle$ extracted by fit, βħω = 5.0")
ylabel!(fitotH, L"$\langle H \rangle$")
xlabel!(fitotH, L"$\eta^2_{max}$")

display(fitotH)
display(fitotx2)

# savefig(fitotH, "..\\imgs_b\\Hfit.png")
# savefig(fitotx2, "..\\imgs_b\\x2fit.png")

# fitx2 = curve_fit(model, esq, vec(df[!,:x2]), df[!,:x2v].^(-2), [0., 0.])
# println("x2")
# println(fitx2.param)
# covx2 = estimate_covar(fitx2)
# println("Errors on params: $(sqrt.(diag(covx2)))")

# fitH = curve_fit(model, esq, vec(df[!,:H]), df[!,:Hv].^(-2), [0., 0.])
# println("H")
# println(fitH.param)
# covH = estimate_covar(fitH)
# println("Errors on params: $(sqrt.(diag(covH)))")
