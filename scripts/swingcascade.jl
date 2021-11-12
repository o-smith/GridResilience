using .NetworkResilience
using PyPlot
using Statistics
using Distributed
using DifferentialEquations
using LinearAlgebra
using PyPlot

q = 0.1
D = 1.0
I = 1.0 
κ = 5.0
αc, fmaxmean, iters = swingcascadebisection(2000, 1e-5, 50, 10, 40, q, 4, κ, I, D)
println(αc)

# println("Starting")

# q = 0.1 
# I = 1.0
# D = 1.0
# κ = 5.0 
# alphas, S, fmaxes, crits = swingcascadevalpha(2000, 50, 25, 25, q, 4, 60, 0.1, 2.5, I, D, κ)

# #Save the values of alpha vs S
# fname = "data/pixelruns/25_25_0_q1_kap5.txt"
# open(fname, "w") do f 
# 	for i=1:length(alphas) 
# 		x = alphas[i] 
# 		y = S[i] 
# 		write(f, "$x $y \n")
# 	end 
# end

# #Save the values of crits
# fname = "data/pixelruns/25_25_0_q1_kap5_hist.txt"
# open(fname, "w") do f 
# 	for i=1:length(crits) 
# 		x = crits[i] 
# 		write(f, "$x \n")
# 	end 
# end 
# q = 0.6
# I = 1.0
# D = 5.0
# κ = 5.0 
# alphas, S, fmaxes, crits = swingcascadevalpha(2000, 50, 25, 25, q, 4, 60, 0.1, 2.5, I, D, κ)

# #Save the values of alpha vs S
# fname = "data/pixelruns/25_25_0_q6_kap5_D5.txt"
# open(fname, "w") do f 
# 	for i=1:length(alphas) 
# 		x = alphas[i] 
# 		y = S[i] 
# 		write(f, "$x $y \n")
# 	end 
# end

# #Save the values of crits
# fname = "data/pixelruns/25_25_0_q6_kap5_D5_hist.txt"
# open(fname, "w") do f 
# 	for i=1:length(crits) 
# 		x = crits[i] 
# 		write(f, "$x \n")
# 	end 
# end 

#hist(crits)
#show() 


# n = 100
# ns = 50 
# nd = 50 
# q = 1.0
# k = 4 
# alpha = 1.0

# #Set up network 
# net, E = smallworldnetwork(n, k, q) 
# sources, sinks = sourcesinklocs(ns, nd, n)
# P = sourcesinkvector(sources, sinks, n, 1.0)
# m = numedges(net)
# nodeset = collect(1:n)
# edgeset = collect(1:m) 
# A = adjacencymatrix(E)
# ψ = rand(2*n) 
# I = 1.0 
# D = 1.0 
# κ = 5.0

# #Time step to find steady state 
# f(u, p, t) = fswing(u, A, P, I, D, κ)
# tspan = (0.0, 100.0) 
# prob = ODEProblem(f, ψ, tspan) 
# sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)

# #Check for convergence 
# ψ = sol[:,end] 
# ω = ψ[1:n]
# θ = ψ[n+1:end]
# x = fsteadystate(θ, A, P, κ) 
# if norm(x,2) > 1e-3
# 	println("Warning: no steady state found.") 
# end

# #Compute the flow and normalise 
# flow_swing = edgepower(θ, E, κ) 
# flowmax_swing = maximum(abs.(flow_swing))
# flow_swing = flow_swing./flowmax_swing 

# #Activate all the edges 
# for i=1:length(net.edgelist)
# 	net.edgelist[i].active = true 
# end

# #Knock out most loaded edge 
# d = argmax(abs.(flow_swing)) 
# edgeset1 = filter(x->x≠d, edgeset) 
# E1 = E[1:end, 1:end .≠ d]
# net.edgelist[d].active = false 

# #Check to see if its broken 
# Adj = adjacencymatrix(E1)
# components, table = connectedcomponents(Adj)

# totedgessurviving = 0
# for i=1:components
# 	println(i)
# 	chunk = table[i] 
# 	n_subset = length(chunk)
# 	nodesubset = zeros(Int64, length(chunk))
# 	ω_subset = zeros(Float64, length(chunk))
# 	θ_subset = zeros(Float64, length(chunk))
# 	ψ_subset = zeros(Float64, 2*n_subset)
# 	for j=1:length(chunk)
# 		nodesubset[j] = nodeset[chunk[j]]
# 		ω_subset[j] = ω[chunk[j]]
# 		θ_subset[j] = θ[chunk[j]]
# 	end
# 	ψ_subset[1:n_subset] = ω_subset
# 	ψ_subset[n_subset+1:end] = θ_subset
# 	edgessurviving = swingfracture!(net, ψ_subset, nodesubset, P, 2.0, alpha, κ, flowmax_swing)
# 	totedgessurviving += edgessurviving
# end

# println(totedgessurviving/m)













