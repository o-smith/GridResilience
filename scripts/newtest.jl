using .NetworkResilience
using DifferentialEquations

# #Test timestepping 
# net, E = smallworldnetwork(16, 4, 0.0) 
# ψ = rand(32).*0.0
# A = adjacencymatrix(E) 
# sources, sinks = sourcesinklocs(1, 15, 16)
# P = sourcesinkvector(sources, sinks, 16, 1.0) 
# obj1 = NetStore(A, P, 5.0) 

# #Find initial steady state
# f(u, p, t) = fswing(u, obj1.adjmat, P, 1.0, 1.0, 5.0)
# tspan = (0.0, 500.0)
# prob = ODEProblem(f, ψ, tspan)  
# sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)
# println(sol[:,end])

# #Test swing cascade bisection 
# ac, fmax, iter = swingcascadebisection(40, 1e-3, 16, 1, 15, 0.0, 4, 5.0, 1.0, 1.0)
# println("ac = ", ac)

n = 60
ns = 30
nd = 30
q = 0.1

κ = 5.0
I = 1.0
D = 1.0 
tol = 5e-4
k = 4
ensemblesize = 500
criticalcouplings = criticalcouplingbisection(ensemblesize, tol, n, ns, nd, q, k, κ, I, D)


fname = "data/critdists/q1_30_30_0.txt"
open(fname, "w") do f 
	for i=1:length(criticalcouplings)
		x = criticalcouplings[i]
		write(f, "$x \n") 
	end
end

