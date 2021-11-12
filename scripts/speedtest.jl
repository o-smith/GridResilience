using .NetworkResilience
using PyPlot

#Make a network 
n = 60
net, E = smallworldnetwork(n, 4, 0.0) 
ψ = rand(2*n) 
A = adjacencymatrix(E) 
sources, sinks = sourcesinklocs(div(n,2), div(n,2), n)
P = sourcesinkvector(sources, sinks, n, 1.0) 
I = 1.0 
D = 1.0 
κ = 5.0

θ = ψ[n+1:end] 
J = mainjac(θ, E, I, κ) 
imshow(J)
colorbar() 
show() 
