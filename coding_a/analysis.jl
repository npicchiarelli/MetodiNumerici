using CSV, DataFrames,Plots

path = "..\\simulations_a\\" 
f1 = path*"clock.csv"

df= CSV.read(f1, DataFrame)

e = df[!,:Energy]
plot(e)