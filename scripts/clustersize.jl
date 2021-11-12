using .NetworkResilience
using PyPlot
using DelimitedFiles

n = 50 
ns = div(n,2) 
nd = div(n,2)
k = 6
ensemblesize = 2500 
s_vec = Array{Float64, 1}([])
q_vec = Array{Float64, 1}([])
d_vec = Array{Float64, 1}([]) 
for q=0:0.01:1
	println(q)
	s = 0.0
	meandegs = 0.0 
	for z=1:ensemblesize
		net, E = smallworldnetwork(n,k,q) 
		meandegs += mean(getdegrees(E))  
		sources, sinks = sourcesinklocs(ns, nd, n)
		A = adjacencymatrix(E) 
		s += clustersize(A, sources) 
	end
	s /= ensemblesize
	meandegs /= ensemblesize 
	push!(s_vec, s)
	push!(q_vec, q)
	push!(d_vec, meandegs) 
end
plot(q_vec, d_vec) #1. .- s_vec./5.9)
show() 

fname = "data/testclusterk6.txt"
open(fname, "w") do f 
	for (i, val) ∈ enumerate(s_vec)
		temp2 = s_vec[i] #1.0 .- s_vec[i]./5.9
		temp = q_vec[i]
		temp3 = d_vec[i]
		write(f, "$temp $temp2 $temp3 \n")
	end
end


# α_vec = Array{Float64, 1}([])
# k_vec = Array{Int64, 1}([])
# for k=2:2:30
# 	println(k)
# 	α, fmax, iter, meandeg = cascadebisection(1, 1e-5, n, ns, nd, q, k, 1.0)
# 	push!(α_vec, α*fmax)
# 	push!(k_vec, k)
# end

# plot(k_vec, α_vec)
# # plot(k_vec, n.*k_vec./(50.0.*(4.0.*k_vec.-4.0)))
# show() 



# net, E = smallworldnetwork(n,k,q) 
# sources, sinks = sourcesinkhalfhalf(ns, nd, n) 
# A = adjacencymatrix(E) 
# s = clustersize(A, sources)
# println(s) 


