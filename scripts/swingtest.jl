using .NetworkResilience
using PyPlot
using DifferentialEquations 
using LinearAlgebra 
using DelimitedFiles 

#To do: make a script to knock out 1 edge, 
#and then make animation of edge power currents 
#whilst also monitoring the deviation norm or something, 
#and plot that as a function of time afterwards,
#and also the residual of fsteadystate.  

#In another script: make script to find how long it takes for the residual 
#to reach zero, or something small. Find this value for different knock out conditions
#and plot against relative load of the knocked out edge. Pruduce a few of these plots for 
#randomly generated lattices of fixed node-type composition.

#Make a network 
n = 30
net, E = smallworldnetwork(n, 2, 0.0) 
ψ = rand(2*n) 
A = adjacencymatrix(E) 
sources, sinks = sourcesinklocs(1, 29, n)
P = sourcesinkvector(sources, sinks, n, 1.0) 

I = 1.0
D = 100.0 
κ = 5.0

f(u, p, t) = fswing(u, A, P, I, D, κ) 

tspan = (0.0, 100.0)
prob = ODEProblem(f, ψ, tspan) 
println("Starting timestepping...")
sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)
ψ = sol[:,end]
θ = sol[n+1:end,end] 

#Knock out most heavily loaded edge 
flow = edgepower(θ, E, κ)
d = argmax(abs.(flow))
E1 = E[1:end, 1:end .≠ d] 
A1 = adjacencymatrix(E1) 
net.edgelist[d].active = false 

#Set up new problem on reduced system 
f(u, p, t) = fswing(u, A1, P, I, D, κ)

#Now start timestepping again 
println("timestepping again...")
tspan = (0.0, 150.0)
prob = ODEProblem(f, ψ, tspan) 
sol = solve(prob,reltol=1e-8,abstol=1e-8) 

# norms = [] 
# anim = @animate for i=1:length(sol.t)
# 	println(i)
# 	y = edgepower(sol[n+1:end,i], E1, κ)
# 	bar(abs.(y))
# 	thetdot = sol[1:n,i]
# 	resid = norm(thetdot, 2)
# 	push!(norms, resid)
# end
# gif(anim, "y_anim_fps15_afterknockout_k4_n30_q1.gif", fps = 15)

thetadots = [] 
thetadoubledots = [] 
resids = [] 
fname = "data/transientresponse/1_29_0_q0_I1_D100.txt"
open(fname, "w") do f 
	for i=1:length(sol.t)
		timenow = sol.t[i]
		ψx = sol[:,i] 
		ωx = ψx[1:n,end]
		θx = ψx[n+1:end] 
		thetadot = norm(ωx,2)
		push!(thetadots, thetadot)  
		ψ_dot = fswing(ψx, A1, P, I, D, κ) 
		θ_doubledot = ψ_dot[1:n]  
		push!(thetadoubledots, norm(θ_doubledot,2))  
		resid = fsteadystate(θx, A1, P, κ) 
		push!(resids, norm(resid,2)) 
		write(f, "$timenow $thetadot \n")
	end 
end

plot(sol.t, thetadots) 
show() 
# savefig("thetadot_kqn_4_1_60")

# plot(sol.t, thetadoubledots) 
# savefig("thetadoubledot_kqn_4_1_60")

# plot(sol.t, resids) 
# savefig("resids_kqn_4_1_60")

# println(sol[:,end])
# ψ_fin = sol[:,end] 
# θ = ψ_fin[n+1:end]
# F = edgepower(θ, E, κ) 
# println(F) 