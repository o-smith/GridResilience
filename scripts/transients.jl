using .NetworkResilience
using PyPlot
using Statistics 
using DifferentialEquations
using LinearAlgebra
using DelimitedFiles


n = 50 
ns = 10
nd = 40 
k = 4 
m = div(n*k,2)
q = 0.0
κ = 5.0
I = 1.0 
D = 10.0 
ensemblesize = 100 

rmat = zeros(m, ensemblesize)
tmat = zeros(m, ensemblesize)
dmat = zeros(m, ensemblesize)

for i=1:ensemblesize
	println(i) 
	r, t, d = centralityvstransient(n, ns, nd, k, q, κ, I, D) 
	rmat[:,i] = r[:]
	tmat[:,i] = t[:]
	dmat[:,i] = d[:]
end

writedlm( "data/trig_vs_cent/10_40_0_q0_I1_D10/rmat.csv",  rmat, ',')
writedlm( "data/trig_vs_cent/10_40_0_q0_I1_D10/tmat.csv",  tmat, ',')
writedlm( "data/trig_vs_cent/10_40_0_q0_I1_D10/dmat.csv",  dmat, ',')

# #Make a network 
# n = 30
# net, E = smallworldnetwork(n, 4, 0.1) 
# ψ = rand(2*n) 
# A = adjacencymatrix(E) 
# sources, sinks = sourcesinklocs(10, 20, n) 
# P = sourcesinkvector(sources, sinks, n, 1.0) 
# I = 1.0 
# D = 1.0 
# κ = 5.0 
# obj1 = CondObj(A, P, κ) 

# #Find the initial steady state 
# f(u, p, t) = fswing(u, obj1.adjmat, P, I, D, κ)  
# tspan = (0.0, 500.0)
# prob = ODEProblem(f, ψ, tspan) 
# println("Starting timestepping...")
# sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)
# ψ0 = sol[:,end]
# θ0 = sol[n+1:end,end] 
# flow = edgepower(θ0, E, κ) 
# meanflow = mean(abs.(flow))

# #Set up callback function 
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
# condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
# affect!(integrator) = terminate!(integrator)

# #Iterate over edges 
# relflow = [] 
# transienttimes = [] 
# distance = [] 
# for d=1:numedges(net) 
# 	println(d) 

# 	#Knock out edge
# 	E1 = E[1:end, 1:end .≠ d] 
# 	A1 = adjacencymatrix(E1) 
# 	obj1.adjmat = A1 
# 	flow_old = flow[1:end .≠ d] 

# 	#Save relative flow on the knocked out edge 
# 	r = abs(flow[d])/meanflow
# 	push!(relflow, r) 

# 	#Set up time stepper and integrate 
# 	u0 = ψ0
# 	tspan2 = (0.0, 5000.0)
# 	prob2 = ODEProblem(f, u0, tspan2) 
# 	cb = DiscreteCallback(condition, affect!) 
# 	sol2 = solve(prob2, reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true,callback=cb)

# 	theta_new = sol2[n+1:end,end]
# 	flow_new = edgepower(theta_new, E1, κ) 
# 	sep = abs.(flow_new .- flow_old) 
# 	push!(distance, norm(sep, 2)) 

# 	#Save length of transient 
# 	push!(transienttimes, sol2.t[end])

# end
# println(transienttimes)
# println(relflow) 

# # plot(relflow[:], transienttimes[:], marker="o", linestyle="")
# plot(relflow[:], distance[:], marker="o", linestyle="", color="k")
# xlabel(" Flow centrality")
# ylabel("Solution distance")
# show() 

# writedlm( "testrmat.csv",  rmat, ',')
# writedlm( "testtmat.csv",  rmat, ',')
# writedlm( "testdmat.csv",  rmat, ',')

# #Object to store information about network for ODE solver 
# mutable struct CondObj
# 	adjmat::Matrix{Float64} 
# 	Pvec::Vector{Float64} 
# 	kappa::Float64 
# end

# #Make a network 
# n = 30
# net, E = smallworldnetwork(n, 4, 1.0) 
# ψ = rand(2*n) 
# A = adjacencymatrix(E) 
# sources, sinks = sourcesinklocs(1, 29, n) 
# P = sourcesinkvector(sources, sinks, n, 1.0) 
# I = 1.0 
# D = 1.0 
# κ = 5.0 
# obj1 = CondObj(A, P, κ) 

# #Find the initial steady state 
# f(u, p, t) = fswing(u, obj1.adjmat, P, I, D, κ)  
# tspan = (0.0, 500.0)
# prob = ODEProblem(f, ψ, tspan) 
# println("Starting timestepping...")
# sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)
# ψ0 = sol[:,end]
# θ0 = sol[n+1:end,end] 
# flow = edgepower(θ0, E, κ) 
# meanflow = mean(abs.(flow))

# #Set up callback function 
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
# condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
# affect!(integrator) = terminate!(integrator)

# #Iterate over edges 
# relflow = [] 
# transienttimes = [] 
# distance = [] 
# for d=1:numedges(net) 
# 	println(d) 

# 	#Knock out edge
# 	E1 = E[1:end, 1:end .≠ d] 
# 	A1 = adjacencymatrix(E1) 
# 	obj1.adjmat = A1 
# 	flow_old = flow[1:end .≠ d] 

# 	#Save relative flow on the knocked out edge 
# 	r = abs(flow[d])/meanflow
# 	push!(relflow, r) 

# 	#Set up time stepper and integrate 
# 	u0 = ψ0
# 	tspan2 = (0.0, 5000.0)
# 	prob2 = ODEProblem(f, u0, tspan2) 
# 	cb = DiscreteCallback(condition, affect!) 
# 	sol2 = solve(prob2, reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true,callback=cb)

# 	theta_new = sol2[n+1:end,end]
# 	flow_new = edgepower(theta_new, E1, κ) 
# 	sep = abs.(flow_new .- flow_old) 
# 	push!(distance, norm(sep, 2)) 

# 	#Save length of transient 
# 	push!(transienttimes, sol2.t[end])

# end
# println(transienttimes)
# println(relflow) 

# plot(relflow[:], transienttimes[:], marker="o", linestyle="")
# plot(relflow[:], distance[:], marker="o", linestyle="", color="k")
# xlabel(" Flow centrality")
# ylabel("Transient time")
# show() 


























# #Set up the callback thing for later 
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
# obj1 = CondObj(A, P, κ) 
# condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
# affect!(integrator) = terminate!(integrator) 

# #Iterate over all possible edge knockouts 
# relflow = [] 
# transienttimes = [] 

# d = 1 
# println(d) 

# #Knock out edge 
# # net.edgelist[:].active .= true 
# E1 = E[1:end, 1:end .≠ d] 
# A1 = adjacencymatrix(E1) 
# # net.edgelist[d].active = false 

# #Save relative flow on the knocked out edge 
# r = abs(flow[d])/meanflow
# push!(relflow, r) 

# #Set up new operator on reduced system 
# f2(u, p, t) = fswing(u, A1, P, I, D, κ) 

# #Set up time stepper 
# u0 = ψ0
# tspan2 = (0.0, 20.0)
# prob2 = ODEProblem(f2, u0, tspan2) 

# #Time step
# obj1.adjmat = A1 
# cb = DiscreteCallback(condition, affect!) 
# sol2 = solve(prob2, reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true, callback=cb) 

# #Save length of transient 
# push!(transienttimes, sol2.t[end]) 


# d = 2 
# println(d) 

# #Knock out edge 
# # net.edgelist[:].active .= true 
# E1 = E[1:end, 1:end .≠ d] 
# A1 = adjacencymatrix(E1) 
# # net.edgelist[d].active = false 

# #Save relative flow on the knocked out edge 
# r = abs(flow[d])/meanflow
# push!(relflow, r) 

# #Set up new operator on reduced system 
# f2(u, p, t) = fswing(u, A1, P, I, D, κ) 

# #Set up timestepper 
# u0 = ψ0
# tspan2 = (0.0, 300.0)
# prob2 = ODEProblem(f2, u0, tspan2) 

# #Timestep
# obj1.adjmat = A1 
# cb = DiscreteCallback(condition, affect!) 
# sol2 = solve(prob2, reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true, callback=cb) 

# #Save length of transient 
# push!(transienttimes, sol2.t[end]) 


# d = 3
# println(d) 

# #Knock out edge 
# # net.edgelist[:].active .= true 
# E1 = E[1:end, 1:end .≠ d] 
# A1 = adjacencymatrix(E1) 
# # net.edgelist[d].active = false 

# #Save relative flow on the knocked out edge 
# r = abs(flow[d])/meanflow
# push!(relflow, r) 

# #Set up new operator on reduced system 
# f2(u, p, t) = fswing(u, A1, P, I, D, κ) 

# #Set up timestepper 
# u0 = ψ0
# tspan2 = (0.0, 300.0)
# prob2 = ODEProblem(f2, u0, tspan2) 

# #Timestep
# obj1.adjmat = A1 
# cb = DiscreteCallback(condition, affect!) 
# sol2 = solve(prob2, reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true, callback=cb) 

# #Save length of transient 
# push!(transienttimes, sol2.t[end]) 

# println(transienttimes)
# # plot(transienttimes) 
# # show() 