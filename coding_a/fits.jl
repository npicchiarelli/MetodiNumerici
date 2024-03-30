using CSV, DataFrames, Dates, Plots, Statistics, LsqFit

size_list = [20, 40, 50, 60, 70, 80]
chimax = []
chimax_v = []
betamax = []

th = plot()
fitp = plot()
constp = plot()



for size in size_list
    f = "..\\simulations_a\\data\\dataL=$size.csv"

    df= CSV.read(f, DataFrame)

    x = df[!,:beta]
    y = df[!,:susc]
    sy = df[!,:susc_v]

    push!(chimax, maximum(y))
    push!(betamax,x[argmax(y)])
    push!(chimax_v, sy[argmax(y)])

    plot!(th,x,y, label = false)
    scatter!(th,x,y, yerr = sy, label = "L = $size", markersize = 1, linestyle = :solid)
end



model(x, p) = p[1].*(x.^(p[2]))

p0 = [1., 7/4]
wt = inv.(chimax_v.^2)
fit = curve_fit(model, size_list, chimax, wt, p0)


println(fit.param)
cov = estimate_covar(fit)
println(cov)
println("Errors on params: $(sqrt.(diag(cov)))")

scatter!(fitp,size_list, chimax, yerr = chimax_v)
xx = collect(LinRange(0,80,1000))
plot!(fitp, xx, model(xx, fit.param))


scatter!(constp, size_list, chimax.*size_list.^(-7/4), yerr = chimax_v.*size_list.^(-7/4))
display(th)
display(constp)
display(fitp)
