using .NetworkResilience

n = 50
k = 4 
q = 0.1
ns = 40
nd = 10
Ptot = 1.0 
net, E = smallworldnetwork(n, k, q) 
sources, sinks = sourcesinklocs(ns, nd, n)
P = noisysourcesinkvector(sources, sinks, n, Ptot)
println(sum(P))
