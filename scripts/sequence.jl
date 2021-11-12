using .NetworkResilience
using PyPlot

#Make ring lattice structure 
n = 10
k = 4
alpha = 0.2
net = Graph(Node[], Edge[], 1, 1) 
for i=1:n 
	addnode!(net) 
end 
for i=1:n 
	for j=1:div(k,2)
		addedge!(net, i, mod1((i+j),n))
	end
end 
E = incidencematrix(net)

#Source-sink locations 
sources = [1,2] 
sinks = [3,4,5,6,7,8,9,10]
# push!(sinks,div(n,2))

T, fracedges, flowseq = cascadesequence(net, E, sources, sinks, alpha, 1.0, 1.0)
println(T)
println(fracedges)
println("Flow sequence:")
for flow ∈ flowseq
	println(flow)
end
