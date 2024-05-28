using CSV, DataFrames,LaTeXStrings, Plots, Printf, Statistics
default(fontfamily = "Computer Modern",
background_color = :transparent,
foreground_color = :black,
margin=5Plots.mm
)

function JKTest(m::Array, blocksize::Int, L::Int)
    length = size(m,1) ÷ blocksize 
    V = L*L

    m_cut = m[1:length*blocksize]

    m_b = reshape(m_cut, (blocksize, length))

    m_abs_blocked = mean(abs, m_b, dims = 1)
    m_abs2_blocked = mean(abs2, m_b, dims = 1)

    jk(a) = let s = sum(a); [s-v for v in a] end

    m_abs_j = jk(m_abs_blocked)./(length-1)
    m_abs2_j = jk(m_abs2_blocked)./(length-1)

    susc = V*(m_abs2_j.-m_abs_j.^2)
    
    return mean(susc), stdm(susc, mean(susc), corrected = false)*sqrt(length-1)
end

p = plot(size = (750,400))
xticks!(p, collect(0:1000:10000))
title!("Standard Deviation of JackKnife Energy as a Function of Block Size")
xlabel!("Block Size [samples]")
ylabel!(L"σ_{Energy}")
plot!(p,[5000, 5000], [0.05, 0.38], lw=2, lc=:black, ls =:dash, label = false)

for (L,β) in [(40,0.8650) (80, 0.8750)]
    Nt = 1000000
    path = "..\\simulations_a\\" 
    v = []
    blocksize_v = collect(Int64, 100:100:10000)

    fname = @sprintf "clock_Nt=%.1eL=%ibeta=%.4f.csv" Nt L β
    f1 = joinpath([path, fname])
    df= CSV.read(f1, DataFrame, types = [Float64, ComplexF64])

    for blocksize in blocksize_v
        magn = df[!,:E]

        means, stds = JKTest(magn,blocksize, L)

        push!(v,stds[1])
        println(blocksize)
    end
    plot!(p,blocksize_v, v, lw = 1.5, label =@sprintf "L=%i beta=%.3f" L β)
end

savefig(p, "..\\imgs_a\\jk.png")
display(p)