using CSV, DataFrames, Dates, Plots, Printf, Statistics, LsqFit, LaTeXStrings, LinearAlgebra
default(fontfamily = "Computer Modern",
background_color = :transparent,
foreground_color = :black,
margin=5Plots.mm
)

function fss(L, beta, y, sy, critexp)
    betac = log(1+sqrt(2))
    return (beta .- betac).*L, y.*L^(-critexp), sy.*L^(-critexp)
end

function fitmax(beta, y, sy, w::Int, plots::Bool)
    am = argmax(y)
    maxy = maximum(y)
    betac = log(1+sqrt(2))
    betac = beta[am]

    parabola(x, p) = p[1] .+ p[2].*(x.-betac).^2

    x_f = beta[am-w:am+w]
    y_f = y[am-w:am+w]
    sy_f = sy[am-w:am+w]

    p0 = [maxy, 7/4]
    wt = inv.(sy_f.^2)
    fit = curve_fit(parabola, x_f, y_f, wt, p0)
    cov = estimate_covar(fit)

    prange = (beta[am] - beta[am-1])*(w+1)
    xx = collect(LinRange(beta[am]-prange, beta[am]+prange, 5000))
    bmax = xx[argmax(parabola(xx, fit.param))]
    berr = xx[2]-xx[1]

    if plots
        fp = plot(titlefontsize = 10)
        xlabel!(fp, "β")
        ylabel!(fp, "χ")
        plot!(fp, xx, parabola(xx, fit.param), linecolor = :black)
        scatter!(fp , x_f, y_f, yerr = sy_f, markersize = 4, markercolor = :black)
        return fit.param[1], sqrt(cov[1]), bmax, berr, fp
    end

    return fit.param[1], sqrt(cov[1]), bmax, berr
end

size_list = [40, 50, 60, 70, 80,]
# size_list = [80,]

chimax = []
chimax_v = []
betamax = []
chimax_fit = []
chimax_fit_v = []
betamax_fit = []
betamax_fit_v = []
fitplots = []

fitp = plot()
constp = plot(legend = false)
Nt = 1000000
Ntstr = @sprintf "%.1e" Nt
th = plot(title = "Specific Heat vs β")
xlabel!(th, L"β")
ylabel!(th, L"C")
fssp = plot(title = "Finite Size Scaling for Specific Heat")
xlabel!(fssp, L"(β-β_c)L^{1/ν}")
ylabel!(fssp, L"CL^{-α/ν}")

for (i,size) in enumerate(size_list)
    f = "..\\simulations_a\\data\\dataL=$size"*"Nt=$Ntstr.csv"

    df= CSV.read(f, DataFrame)

    beta = df[!,:beta]
    y = df[!,:susc]
    sy = df[!,:susc_v]
    am = argmax(y)
    maxy = maximum(y)

    fmax, maxerr, bmaxf, berrf, fitplot = fitmax(beta, y, sy, 3, true)
    title!(fitplot, "L=$size",)
    push!(fitplots, fitplot)
    
    beta_fss, chi_fss, schi_fss = fss(size, beta, y, sy, 0)

    println("$size:max = $(maxy), argmax = $(am) @ beta = $(beta[am])")
    println("Max from fit: $fmax ± $maxerr @ $bmaxf")
    push!(chimax, maxy)
    push!(betamax,beta[am])
    push!(chimax_v, sy[am])

    push!(chimax_fit, fmax)
    push!(chimax_fit_v, maxerr)
    push!(betamax_fit, bmaxf)
    push!(betamax_fit_v, 2berrf)

    color = palette(:tab10)[i]
    # plot!(th,beta,y, label = false, linecolor = color)
    scatter!(th,beta,y, yerr = sy, label = "L = $size", markersize = 4, markercolor = color)

    # plot!(fssp,beta_fss,chi_fss, label = false, linecolor = color)
    scatter!(fssp,beta_fss,chi_fss, yerr = schi_fss, label = "L = $size", markersize = 4, markercolor = color)
end

bigplot = plot(fitplots..., layout = (5,1), size=(750, 1000), legend = false, plot_title = "Parabolic fits with 7 points near maximum")

betamodel(x,p) = p[1] .+ p[2].*x.^(-1/p[3])
p0 = [log(1+sqrt(2)), 1., 1.]
wt = fill(0.001, size(betamax))
betafit = curve_fit(betamodel, size_list, betamax, wt.^(-2), p0)

println("Beta")
println(betafit.param)
cov = estimate_covar(betafit)
println("Errors on params: $(sqrt.(diag(cov)))")

model(x, p) = p[1].*(x.^(p[2]))

p0 = [1., 7/4]
wt = inv.(chimax_v.^2)
fit = curve_fit(model, size_list, chimax, wt, p0)

println(fit.param)
cov = estimate_covar(fit)
println("Errors on params: $(sqrt.(diag(cov)))")

scatter!(fitp,size_list, chimax, yerr = chimax_v, label = "Sampled Max")
xx = collect(LinRange(0,80,1000))
plot!(fitp, xx, model(xx, fit.param), label = "Best Fit")
title!(fitp, "Critical exponent fit, Sampled Max")
ylabel!(fitp, L"χ_{max}")
xlabel!(fitp, L"L")

scatter!(constp, size_list, chimax.*size_list.^(-7/4), yerr = chimax_v.*size_list.^(-7/4))
title!(constp, "Finite Size Scaling for Susceptivity Peak")
xlabel!(constp, L"L")
ylabel!(constp, L"χ_{max}L^{-γ/ν}")

fitf = curve_fit(model, size_list, chimax, chimax_fit_v.^(-2), p0)
println(fitf.param)
covf = estimate_covar(fitf)
println(covf)
println("Errors on params: $(sqrt.(diag(covf)))")

fitpf = plot(title = "Critical exponent fit, Fit Max")
scatter!(fitpf,size_list, chimax_fit, yerr = chimax_fit_v, label = "Max from Fit")
plot!(fitpf, xx, model(xx, fitf.param), label = "Best Fit")
ylabel!(fitpf, L"χ_{max}")
xlabel!(fitpf, L"L")

constpf = plot(legend = false)
scatter!(constpf, size_list, chimax_fit.*size_list.^(-7/4), yerr = chimax_fit_v.*size_list.^(-7/4))
title!(constpf, "Finite Size Scaling for Susceptivity Peak, Fit")
xlabel!(constpf, L"L")
ylabel!(constpf, L"χ_{max}L^{-γ/ν}")

crit = plot(title = "Comparison between maximum methods, γ")
scatter!(crit, [1,2], [fit.param[2], fitf.param[2]], yerr = sqrt.([cov[4], covf[4]]), ms = 5, markershape=:square, markercolor=:black, label = false)
plot!(crit, [0.5, 2.5], [7/4, 7/4], lw=2, lc=:black, ls =:dash, label = "Theory") 
xticks!(crit, [1,2], ["Simple Max", "Fit Max"])
ylabel!(crit, "γ")

# display(bigplot)
# display(th)
# display(fssp)

display(constp)
# display(fitp)

display(constpf)
# display(fitpf)
# display(crit)

# savefig(bigplot, "..\\imgs_a\\bigplot.png")
# savefig(th, "..\\imgs_a\\spec.png")
# savefig(fssp, "..\\imgs_a\\spec_fss.png")

# savefig(constp, "..\\imgs_a\\constp.png")
# savefig(fitp, "..\\imgs_a\\fitp.png")

# savefig(constpf, "..\\imgs_a\\constpf.png")
# savefig(fitpf, "..\\imgs_a\\fitpf.png")
# savefig(crit, "..\\imgs_a\\crit.png")
