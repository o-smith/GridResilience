using .NetworkResilience
using PyPlot
using Statistics 
using DifferentialEquations
using LinearAlgebra

#Make test network 
#Make a network 
n = 30
net, E = smallworldnetwork(n, 4, 1.0) 
ψ = rand(2*n) 
A = adjacencymatrix(E) 
sources = [1]
sinks = [16] 
P = sourcesinkvector(sources, sinks, n, 1.0) 
I = 1.0 
D = 1.0 
κ = 5.0

#Solve initial condition 
f(u, p, t) = fswing(u, A, P, I, D, κ)
tspan = (0.0, 100.0) 
prob = ODEProblem(f, ψ, tspan) 
println("Starting timestepping...")
sol = solve(prob,reltol=1e-8,abstol=1e-8) #,save_everystep=false,dense=false,save_end=true)
ψ = sol[:,end]
θ = sol[n+1:end,end] 
flow1 = edgepower(θ, E, κ) 

m = numedges(net)
X = zeros(Float64, m) .+ 1.0

L = weightedlaplacianmatrix(E, X) 
theta = nodevoltages(L, P)
F = edgecurrents(theta, E, X)

println(abs.(flow1[:]))
println(abs.(F[:]))

plot(abs.(flow1[:]), color="k", linewidth=3.0)
plot(abs.(F), color="r", linestyle="--", linewidth=2.0)
show() 

# #Knock out an edge 
# d = 10 
# E1 = E[1:end, 1:end .≠ d] 
# A1 = adjacencymatrix(E1) 

# #Set up new problem on reduced system 
# f(u, p, t) = fswing(u, A1, P, I, D, κ)

# #Now start timestepping again 
# println("timestepping again...")
# tspan = (0.0, 300.0)
# prob = ODEProblem(f, ψ, tspan) 
# sol = solve(prob,reltol=1e-8,abstol=1e-8) 

# resids = [] 
# for i=1:length(sol.t)
# 	ψx = sol[:,i] 
# 	θx = ψx[n+1:end] 
# 	resid = fsteadystate(θx, A1, P, κ) 
# 	push!(resids, norm(resid,2)) 
# end 

# plot(sol.t, resids) 
# show() 

# #Now do it again 
# mutable struct CondObj
# 	adjmat::Matrix{Float64} 
# 	Pvec::Vector{Float64} 
# 	kappa::Float64 
# end
# function conditiontemp(u1, t1, integrator, ob::CondObj)   
# 	n = div(length(u1),2) 
# 	θ = u1[n+1:end,end]
# 	resid = fsteadystate(θ, ob.adjmat, ob.Pvec, ob.kappa)
# 	if norm(resid,2) < 0.001
# 		return true 
# 	else
# 		return false
# 	end
# end
# obj1 = CondObj(A1, P, κ) 
# condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
# affect!(integrator) = terminate!(integrator)
# cb = DiscreteCallback(condition, affect!) 


# f(u, p, t) = fswing(u, A1, P, I, D, κ)
# tspan = (0.0, 300.0)
# prob = ODEProblem(f, ψ, tspan) 
# sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true, callback=cb)
# println(sol.t[end])












