using .NetworkResilience
using PyPlot

#System parameters
alpha = 5.0 
q = 0.1
n = 50
k = 4
ns = 10
nd = 10
Ptot = 20000.0 

#Make network
net, E0 = smallworldnetwork(n, k, q) 
sources, sinks = sourcesinklocs(ns, nd, n)
P = sourcesinkvector(sources, sinks, n, Ptot) 
n = numnodes(net) 
m = numedges(net)
nodeset = collect(1:n)
edgeset = collect(1:m)  
X = zeros(Float64, m) 
X .+= 0.5 

#Solve the DC flow equations
L = weightedlaplacianmatrix(E0, 1.0./X)  
theta = nodevoltages(L, P)
F = edgecurrents(theta, E0, X)

#Knock out most loaded edge 
d = argmax(abs.(F)) 
edgeset1 = filter(x->x≠d, edgeset) 
E1 = E0[1:end, 1:end .≠ d]
net.edgelist[d].active = false 

#Check to see if its broken 
Adj = adjacencymatrix(E1)
components, table = connectedcomponents(Adj) 

info = newrecorder(m, ns, Ptot, nd, sum(abs.(F.*X)))

# #Stuff needed for function call 
# info = [0.0 Float64(length(edgeset))] 
simtime = 0

for i=1:components
	chunk = table[i] 
	nodesubset = zeros(Int64, length(chunk))
	for j=1:length(chunk)
		nodesubset[j] = nodeset[chunk[j]]
	end
	global info = fracture!(net, table[i], P, X, F, info, simtime, alpha)
end

v = hcat(info.simulationtime, info.liveedges, info.weightedlivesources)
println(v)




