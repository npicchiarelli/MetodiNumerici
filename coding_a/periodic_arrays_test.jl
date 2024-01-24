M = [[1 2 3]; [4 5 6]; [7 8 9]]
ci = CartesianIndices(M)
I1 = oneunit(ci[1])
for i in ci 
    @show i_v = vcat(getindex(i,1) ,getindex(i,2))
end

# lrud = [CartesianIndex(0,1) CartesianIndex(0,-1) CartesianIndex(1,0) CartesianIndex(-1,0)]

# i = CartesianIndex(1,1)
# j = CartesianIndex(2,3)
# a = [i j]

# ind = findall(<(5), M)

# f = hcat([i[1] for i in intuple], [i[2] for i in ind])
# println(f)
# replace!(f, 1=>5, 2=>4)
# print(f)
# ff = CartesianIndex.(f[:,1], f[:,2])
# print(ff)

# lrud = [[0 1]; [0 -1]; [1 0]; [-1 0]] 
# i_v = hcat(getindex(i,1) ,getindex(i,2))
# nn = i_v .+ lrud
# print(nn)