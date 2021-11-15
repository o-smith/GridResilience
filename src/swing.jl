function sinematrix(s::Vector{Float64}, c::Vector{Float64})
	return s.*c' .- c.*s' 
end


function cosmatrix(s::Vector{Float64}, c::Vector{Float64})
	return c.*c' + s.*s' 
end


function cmatrix(κ::Float64, A::Matrix{Float64}, S::Matrix{Float64})
	return κ.*A.*S 
end


function omegadot(ω::Vector{Float64}, θ::Vector{Float64}, A::Matrix{Float64}, P::Vector{Float64}, I::Float64, D::Float64, κ::Float64)
	s = sin.(θ)
	c = cos.(θ)
	smat = sinematrix(s, c) 
	C = cmatrix(κ, A, smat) 
	return (P .- D.*ω .- C*ones(length(ω)))./I 
end


function fswing(ψ::Vector{Float64}, A::Matrix{Float64}, P::Vector{Float64}, I::Float64, D::Float64, κ::Float64)
	m = length(ψ)
	n = div(m,2) 
	ω = ψ[1:n]
	θ = ψ[n+1:end]
	ω_dot = omegadot(ω, θ, A, P, I, D, κ) 
	return vcat(ω_dot, ω) 
end


function fsteadystate(θ::Vector{Float64}, A::Matrix{Float64}, P::Vector{Float64}, κ::Float64)
	n = length(θ)   
	s = sin.(θ)
	c = cos.(θ)
	smat = sinematrix(s, c) 
	C = cmatrix(κ, A, smat) 
	return P .- C*ones(n) 
end 


function mainjac(θ::Vector{Float64}, E::Matrix{Float64}, I::Float64, κ::Float64)
	weights = κ*cos.(E'*θ) 
	W = Diagonal(weights)
	return -E*W*E'  
end


function dfswing(ψ::Vector{Float64}, A::Matrix{Float64}, I::Float64, D::Float64, κ::Float64)
	m = length(ψ) 
	n = div(m,2) 
	θ = ψ[n+1:end] 
	DF = zeros(m,m)
	DF[1:n,1:n] = D.*Matrix{Float64}(I, n, n)./I
	DF[1:n,n+1:end] = mainjac(θ, A, I, κ)./I 
	DF[n+1:end,1:n] = Matrix{Float64}(I, n, n)
	return DF 
end


function edgepower(θ::Vector{Float64}, E::Matrix{Float64}, κ::Float64)
	Δθ = transpose(E)*θ
	return κ.*sin.(Δθ) 
end


function centralityvstransient(n::Int64, ns::Int64, nd::Int64, k::Int64, q::Float64, κ::Float64, I::Float64, D::Float64) 

	#Set up network
	net, E = smallworldnetwork(n, k, q) 
	ψ = rand(2*n)
	A = adjacencymatrix(E) 
	sources, sinks = sourcesinklocs(ns, nd, n)
	P = sourcesinkvector(sources, sinks, n, 1.0) 
	obj1 = NetStore(A, P, κ) 

	#Find initial steady state
	f(u, p, t) = fswing(u, obj1.adjmat, P, I, D, κ)
	tspan = (0.0, 500.0)
	prob = ODEProblem(f, ψ, tspan)  
	sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)
	ψ0 = sol[:,end]
	θ0 = sol[n+1:end,end] 
	flow = edgepower(θ0, E, κ) 
	meanflow = mean(abs.(flow))

	#Set up callback function 
	function conditiontemp(u1, t1, integrator, ob::NetStore)   
		n = div(length(u1),2) 
		θ = u1[n+1:end,end]
		resid = fsteadystate(θ, ob.adjmat, ob.pvec, ob.kappa)
		if norm(resid,2) < 0.001
			return true 
		else
			return false
		end
	end
	condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
	affect!(integrator) = terminate!(integrator)

	m = numedges(net)
	relflow = zeros(m) 
	transienttimes = zeros(m) 
	distances = zeros(m) 
	for d=1:m 

		#Delete edge and find its load 
		E1 = E[1:end, 1:end .≠ d] 
		A1 = adjacencymatrix(E1) 
		obj1.adjmat = A1 
		relflow[d] = abs(flow[d])/meanflow
		flow_old = flow[1:end .≠ d] 

		#Time step on reduced network 
		u0 = ψ0
		tspan2 = (0.0, 5000.0)
		prob2 = ODEProblem(f, u0, tspan2) 
		cb = DiscreteCallback(condition, affect!) 
		sol2 = solve(prob2, reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true,callback=cb)
		
		#Find new flow and measure its distance 
		theta_new = sol2[n+1:end,end]
		flow_new = edgepower(theta_new, E1, κ) 
		sep = abs.(flow_new .- flow_old) 
		transienttimes[d] = sol2.t[end]
		distances[d] = norm(sep, 2) 
	end
	return relflow, transienttimes, distances  

end 


function swingfracture!(net::Graph, ψ::Vector{Float64}, nodeset::Vector{Int64}, P::Vector{Float64}, synctol::Float64, alpha::Float64, kappa::Float64, maxflow::Float64)

	tol = 1e-5

	#Check the status of P, see if has any sources or sinks. If not, halt recursion
	nodecardinality = length(nodeset)
	P1 = zeros(Float64, nodecardinality) 
	sourcecounter = 0 
	sinkcounter = 0 
	for i=1:nodecardinality
		vi = nodeset[i] 
		P1[i] = P[vi] 
		if P1[i] > tol 
			sourcecounter += 1 
		end
		if P1[i] < -tol 
			sinkcounter += 1
		end
	end
	if (sourcecounter == 0) || (sinkcounter == 0)
		return 0 
	end

	#Balance P1 homogeneously 
	delta = sum(P1)/2.0 
	individualdelta_sink = delta/Float64(sinkcounter)
	individualdelta_source = delta/Float64(sourcecounter)
	for i=1:nodecardinality
		if P1[i] < -tol 
			P1[i] -= individualdelta_sink  
		end 
		if P1[i] > tol
			P1[i] -= individualdelta_source
		end
	end 

	# #Balance P1 - if there is power excess, increase loads
	# delta = sum(P1) 
	# individualdelta = delta/Float64(sinkcounter)
	# for i=1:nodecardinality
	# 	if P1[i] < -tol 
	# 		P1[i] -= individualdelta
	# 	end 
	# end  

	#Find which edges are connected to the node-set
	edgeset = Array{Int64, 1}([]) 
	for vi ∈ nodeset
		v = net.nodelist[vi]
		for e ∈ v.outedges
			if (e ∉ edgeset) && (net.edgelist[e].active == true)
				push!(edgeset, e)
			end 	
		end
		for e ∈ v.inedges 
			if (e ∉ edgeset) && (net.edgelist[e].active == true)
				push!(edgeset, e)
			end 	
		end 
	end
	edgecardinality = length(edgeset)

	#Construct E1 
	E1 = zeros(Float64, nodecardinality, edgecardinality) 
	for i=1:nodecardinality
		vi = nodeset[i] 
		v = net.nodelist[vi]
		for j=1:edgecardinality
			if edgeset[j] ∈ v.outedges
				E1[i,j] = 1.0
			end
			if edgeset[j] ∈ v.inedges
				E1[i,j] = -1.0
			end 
		end		
	end

	#Set up the storage object 
	A1 = adjacencymatrix(E1) 
	obj1 = NetStore2(A1, P1, E1, kappa, synctol, alpha, maxflow, true, false) 

	#Set up callback function 
	function conditiontemp(u1, t1, integrator, ob::NetStore2)   
		n = div(length(u1),2) 
		ω = u1[1:n]
		θ = u1[n+1:end]

		if (t1 > 0.0)
			#Has network synchronised? 
			if norm(ω,2) > ob.synctol 
				ob.sync = false
				return true 
			end 

			#Have any edges tripped the capacity 
			flow = edgepower(θ, ob.incidence, ob.kappa)./ob.maxflow
			for i=1:length(flow)
				if abs(flow[i]) > ob.alpha
					ob.trip = true 
					return true 
				end
			end

			#Has the solution converged
			resid = fsteadystate(θ, ob.adjmat, ob.pvec, kappa)
			if norm(resid,2) < 1e-6
				return true 
			end 
		end 

		#If the network is still syncing and no edge trips, continue 
		return false 
	end
	condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
	affect!(integrator) = terminate!(integrator)

	#Time step 
	f(u, p, t) = fswing(u, obj1.adjmat, P1, 1.0, 1.0, kappa)
	tspan = (0.0, 500.0)
	prob = ODEProblem(f, ψ, tspan) 
	cb = DiscreteCallback(condition, affect!) 
	sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true,callback=cb)

	#Get the state vector, split into theta and omega 
	ψ_1 = sol[:,end]
	ω_1 = ψ_1[1:nodecardinality]
	θ_1 = ψ_1[nodecardinality+1:end]

	#If the solution has not converged, then return
	if obj1.sync == false 
		return 0 
	end 

	#If the solution has converged with no further overloads, then return 
	if (obj1.trip == false) && (obj1.sync == true)  
		return edgecardinality
	end 

	#If there are overloads, but omega still in acceptable range, 
	#then compute number of new connected components and continue 
	#recursion on each of them
	#First construct new edgeset and incidence matrix 
	overloadcounter = 0
	survivors = Array{Int64, 1}([]) 
	flow = edgepower(θ_1, E1, kappa)./maxflow
	for i=1:edgecardinality
		ei = edgeset[i]
		if abs(flow[i]) > alpha #*abs(Finit[ei]) 
			net.edgelist[ei].active = false
			overloadcounter += 1 
		else
			push!(survivors, i)
		end 
	end
	newedgecard = length(survivors) 
	E2 = zeros(Float64, nodecardinality, newedgecard) 
	for j=1:newedgecard
		E2[:,j] = E1[:,survivors[j]] 
	end

	#Check to see if its broken 
	Adj = adjacencymatrix(E2)
	components, table = connectedcomponents(Adj)

	#Clear any superfluous data
	E1 = nothing 
	E2 = nothing
	Adj = nothing
	obj1 = nothing 

	#Call fracture function on each of the components, passing down info
	descendent_edges = 0
	for i=1:components
		chunk = table[i] 
		n_subset = length(chunk)
		nodesubset = zeros(Int64, length(chunk))
		ω_subset = zeros(Float64, length(chunk))
		θ_subset = zeros(Float64, length(chunk))
		ψ_subset = zeros(Float64, 2*n_subset)
		for j=1:length(chunk)
			nodesubset[j] = nodeset[chunk[j]] 
			ω_subset[j] = ω_1[chunk[j]]
			θ_subset[j] = θ_1[chunk[j]]
		end
		ψ_subset[1:n_subset] = ω_subset
		ψ_subset[n_subset+1:end] = θ_subset
		descendent_edges += swingfracture!(net, ψ_subset, nodesubset, P, synctol, alpha, kappa, maxflow) 
	end 

	return descendent_edges
end


function swingfracturewithtime!(net::Graph, ψ::Vector{Float64}, nodeset::Vector{Int64}, P::Vector{Float64}, synctol::Float64, alpha::Float64, kappa::Float64, maxflow::Float64, D::Float64)

	tol = 1e-5

	#Check the status of P, see if has any sources or sinks. If not, halt recursion
	nodecardinality = length(nodeset)
	P1 = zeros(Float64, nodecardinality) 
	sourcecounter = 0 
	sinkcounter = 0 
	for i=1:nodecardinality
		vi = nodeset[i] 
		P1[i] = P[vi] 
		if P1[i] > tol 
			sourcecounter += 1 
		end
		if P1[i] < -tol 
			sinkcounter += 1
		end
	end
	if (sourcecounter == 0) || (sinkcounter == 0)
		return 0, 0.0 
	end

	#Balance P1 homogeneously 
	delta = sum(P1)/2.0 
	individualdelta_sink = delta/Float64(sinkcounter)
	individualdelta_source = delta/Float64(sourcecounter)
	for i=1:nodecardinality
		if P1[i] < -tol 
			P1[i] -= individualdelta_sink  
		end 
		if P1[i] > tol
			P1[i] -= individualdelta_source
		end
	end 

	# #Balance P1 - if there is power excess, increase loads
	# delta = sum(P1) 
	# individualdelta = delta/Float64(sinkcounter)
	# for i=1:nodecardinality
	# 	if P1[i] < -tol 
	# 		P1[i] -= individualdelta
	# 	end 
	# end  

	#Find which edges are connected to the node-set
	edgeset = Array{Int64, 1}([]) 
	for vi ∈ nodeset
		v = net.nodelist[vi]
		for e ∈ v.outedges
			if (e ∉ edgeset) && (net.edgelist[e].active == true)
				push!(edgeset, e)
			end 	
		end
		for e ∈ v.inedges 
			if (e ∉ edgeset) && (net.edgelist[e].active == true)
				push!(edgeset, e)
			end 	
		end 
	end
	edgecardinality = length(edgeset)

	#Construct E1 
	E1 = zeros(Float64, nodecardinality, edgecardinality) 
	for i=1:nodecardinality
		vi = nodeset[i] 
		v = net.nodelist[vi]
		for j=1:edgecardinality
			if edgeset[j] ∈ v.outedges
				E1[i,j] = 1.0
			end
			if edgeset[j] ∈ v.inedges
				E1[i,j] = -1.0
			end 
		end		
	end

	#Set up the storage object 
	A1 = adjacencymatrix(E1) 
	obj1 = NetStore2(A1, P1, E1, kappa, synctol, alpha, maxflow, true, false) 

	#Set up callback function 
	function conditiontemp(u1, t1, integrator, ob::NetStore2)   
		n = div(length(u1),2) 
		ω = u1[1:n]
		θ = u1[n+1:end]

		if (t1 > 0.0)
			#Has network synchronised? 
			if norm(ω,2) > ob.synctol 
				ob.sync = false
				return true 
			end 

			#Have any edges tripped the capacity 
			flow = edgepower(θ, ob.incidence, ob.kappa)./ob.maxflow
			for i=1:length(flow)
				if abs(flow[i]) > ob.alpha
					ob.trip = true 
					return true 
				end
			end

			#Has the solution converged
			resid = fsteadystate(θ, ob.adjmat, ob.pvec, kappa)
			if norm(resid,2) < 1e-6
				return true 
			end 
		end 

		#If the network is still syncing and no edge trips, continue 
		return false 
	end
	condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
	affect!(integrator) = terminate!(integrator)

	#Time step 
	f(u, p, t) = fswing(u, obj1.adjmat, P1, 1.0, D, kappa)
	tspan = (0.0, 500.0)
	prob = ODEProblem(f, ψ, tspan) 
	cb = DiscreteCallback(condition, affect!) 
	sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true,callback=cb)

	#Get the state vector, split into theta and omega 
	ψ_1 = sol[:,end]
	ω_1 = ψ_1[1:nodecardinality]
	θ_1 = ψ_1[nodecardinality+1:end]

	#If the solution has not converged, then return
	if obj1.sync == false 
		return 0, 0.0 
	end 

	#If the solution has converged with no further overloads, then return 
	if (obj1.trip == false) && (obj1.sync == true)  
		return edgecardinality, sol.t[end] 
	end 

	#If there are overloads, but omega still in acceptable range, 
	#then compute number of new connected components and continue 
	#recursion on each of them
	#First construct new edgeset and incidence matrix 
	overloadcounter = 0
	survivors = Array{Int64, 1}([]) 
	flow = edgepower(θ_1, E1, kappa)./maxflow
	for i=1:edgecardinality
		ei = edgeset[i]
		if abs(flow[i]) > alpha #*abs(Finit[ei]) 
			net.edgelist[ei].active = false
			overloadcounter += 1 
		else
			push!(survivors, i)
		end 
	end
	newedgecard = length(survivors) 
	E2 = zeros(Float64, nodecardinality, newedgecard) 
	for j=1:newedgecard
		E2[:,j] = E1[:,survivors[j]] 
	end

	#Check to see if its broken 
	Adj = adjacencymatrix(E2)
	components, table = connectedcomponents(Adj)

	#Clear any superfluous data
	E1 = nothing 
	E2 = nothing
	Adj = nothing
	obj1 = nothing 

	#Call fracture function on each of the components, passing down info
	descendent_edges = 0
	times = []
	for i=1:components
		chunk = table[i] 
		n_subset = length(chunk)
		nodesubset = zeros(Int64, length(chunk))
		ω_subset = zeros(Float64, length(chunk))
		θ_subset = zeros(Float64, length(chunk))
		ψ_subset = zeros(Float64, 2*n_subset)
		for j=1:length(chunk)
			nodesubset[j] = nodeset[chunk[j]] 
			ω_subset[j] = ω_1[chunk[j]]
			θ_subset[j] = θ_1[chunk[j]]
		end
		ψ_subset[1:n_subset] = ω_subset
		ψ_subset[n_subset+1:end] = θ_subset
		xz, tz = swingfracturewithtime!(net, ψ_subset, nodesubset, P, synctol, alpha, kappa, maxflow, D) 
		descendent_edges += xz
		push!(times, tz)
	end 
	tprime = sol.t[end] + maximum(times) 
	return descendent_edges, tprime
end


function swingfracturewithbreakdown!(net::Graph, ψ::Vector{Float64}, nodeset::Vector{Int64}, P::Vector{Float64}, synctol::Float64, alpha::Float64, kappa::Float64, maxflow::Float64)

	tol = 1e-5

	#Check the status of P, see if has any sources or sinks. If not, halt recursion
	nodecardinality = length(nodeset)
	P1 = zeros(Float64, nodecardinality) 
	sourcecounter = 0 
	sinkcounter = 0 
	for i=1:nodecardinality
		vi = nodeset[i] 
		P1[i] = P[vi] 
		if P1[i] > tol 
			sourcecounter += 1 
		end
		if P1[i] < -tol 
			sinkcounter += 1
		end
	end

	#Find which edges are connected to the node-set
	edgeset = Array{Int64, 1}([]) 
	for vi ∈ nodeset
		v = net.nodelist[vi]
		for e ∈ v.outedges
			if (e ∉ edgeset) && (net.edgelist[e].active == true)
				push!(edgeset, e)
			end 	
		end
		for e ∈ v.inedges 
			if (e ∉ edgeset) && (net.edgelist[e].active == true)
				push!(edgeset, e)
			end 	
		end 
	end
	edgecardinality = length(edgeset)
	if (sourcecounter == 0) || (sinkcounter == 0)
		return 0, 0, 0, edgecardinality 
	end

	#Balance P1 homogeneously 
	delta = sum(P1)/2.0 
	individualdelta_sink = delta/Float64(sinkcounter)
	individualdelta_source = delta/Float64(sourcecounter)
	for i=1:nodecardinality
		if P1[i] < -tol 
			P1[i] -= individualdelta_sink  
		end 
		if P1[i] > tol
			P1[i] -= individualdelta_source
		end
	end 


	#Construct E1 
	E1 = zeros(Float64, nodecardinality, edgecardinality) 
	for i=1:nodecardinality
		vi = nodeset[i] 
		v = net.nodelist[vi]
		for j=1:edgecardinality
			if edgeset[j] ∈ v.outedges
				E1[i,j] = 1.0
			end
			if edgeset[j] ∈ v.inedges
				E1[i,j] = -1.0
			end 
		end		
	end

	#Set up the storage object 
	A1 = adjacencymatrix(E1) 
	obj1 = NetStore2(A1, P1, E1, kappa, synctol, alpha, maxflow, true, false) 

	#Set up callback function 
	function conditiontemp(u1, t1, integrator, ob::NetStore2)   
		n = div(length(u1),2) 
		ω = u1[1:n]
		θ = u1[n+1:end]

		if (t1 > 0.0)
			#Has network synchronised? 
			if norm(ω,2) > ob.synctol 
				ob.sync = false
				return true 
			end 

			#Have any edges tripped the capacity 
			flow = edgepower(θ, ob.incidence, ob.kappa)./ob.maxflow
			for i=1:length(flow)
				if abs(flow[i]) > ob.alpha
					ob.trip = true 
					return true 
				end
			end

			#Has the solution converged
			resid = fsteadystate(θ, ob.adjmat, ob.pvec, kappa)
			if norm(resid,2) < 1e-6
				return true 
			end 
		end 

		#If the network is still syncing and no edge trips, continue 
		return false 
	end
	condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
	affect!(integrator) = terminate!(integrator)

	#Time step 
	f(u, p, t) = fswing(u, obj1.adjmat, P1, 1.0, 1.0, kappa)
	tspan = (0.0, 500.0)
	prob = ODEProblem(f, ψ, tspan) 
	cb = DiscreteCallback(condition, affect!) 
	sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true,callback=cb)

	#Get the state vector, split into theta and omega 
	ψ_1 = sol[:,end]
	ω_1 = ψ_1[1:nodecardinality]
	θ_1 = ψ_1[nodecardinality+1:end]

	#If the solution has not converged, then return
	if obj1.sync == false 
		return 0, 0, edgecardinality, 0
	end 

	#If the solution has converged with no further overloads, then return 
	if (obj1.trip == false) && (obj1.sync == true)  
		return edgecardinality, 0, 0, 0 
	end 

	#If there are overloads, but omega still in acceptable range, 
	#then compute number of new connected components and continue 
	#recursion on each of them
	#First construct new edgeset and incidence matrix 
	overloadcounter = 0
	survivors = Array{Int64, 1}([]) 
	flow = edgepower(θ_1, E1, kappa)./maxflow
	for i=1:edgecardinality
		ei = edgeset[i]
		if abs(flow[i]) > alpha 
			net.edgelist[ei].active = false
			overloadcounter += 1 
		else
			push!(survivors, i)
		end 
	end
	newedgecard = length(survivors) 
	E2 = zeros(Float64, nodecardinality, newedgecard) 
	for j=1:newedgecard
		E2[:,j] = E1[:,survivors[j]] 
	end

	#Check to see if its broken 
	Adj = adjacencymatrix(E2)
	components, table = connectedcomponents(Adj)

	#Clear any superfluous data
	E1 = nothing 
	E2 = nothing
	Adj = nothing
	obj1 = nothing 

	#Call fracture function on each of the components, passing down info
	descendent_edges = 0
	descendent_overloads = 0
	descendent_desyncs = 0 
	descendent_unbalances = 0 
	for i=1:components
		chunk = table[i] 
		n_subset = length(chunk)
		nodesubset = zeros(Int64, length(chunk))
		ω_subset = zeros(Float64, length(chunk))
		θ_subset = zeros(Float64, length(chunk))
		ψ_subset = zeros(Float64, 2*n_subset)
		for j=1:length(chunk)
			nodesubset[j] = nodeset[chunk[j]] 
			ω_subset[j] = ω_1[chunk[j]]
			θ_subset[j] = θ_1[chunk[j]]
		end
		ψ_subset[1:n_subset] = ω_subset
		ψ_subset[n_subset+1:end] = θ_subset
		xz, oz, dz, uz = swingfracturewithbreakdown!(net, ψ_subset, nodesubset, P, synctol, alpha, kappa, maxflow) 
		descendent_edges += xz
		descendent_overloads += oz 
		descendent_desyncs += dz 
		descendent_unbalances += uz 
	end 
	return descendent_edges, descendent_overloads, descendent_desyncs, descendent_unbalances 
end





function swingcascadeandtvalpha(ensemblesize::Int64, n::Int64, ns::Int64, nd::Int64, q::Float64, 
								k::Int64, alphares::Int64, alphamin::Float64, alphamax::Float64,
								I::Float64, D::Float64, κ::Float64)

	alphas = collect(range(alphamin, length=alphares, stop=alphamax))
	edgemat = zeros(Float64, ensemblesize, alphares) 
	tmat = zeros(Float64, ensemblesize, alphares)
	fmaxes = zeros(Float64, ensemblesize)
	critvalues = zeros(Float64, ensemblesize) 

	println("Starting")
	for z=1:ensemblesize
		println("Ensemble = ", z)

		net, E = smallworldnetwork(n, k, q)
		# net, E = scalefreenetwork(n, 2*n, 3) 
		sources, sinks = sourcesinklocs(ns, nd, n)
		P = sourcesinkvector(sources, sinks, n, 1.0)
		m = numedges(net)
		nodeset = collect(1:n)
		edgeset = collect(1:m) 
		A = adjacencymatrix(E)
		ψ = rand(2*n) 

		#Time step to find steady state 
		f(u, p, t) = fswing(u, A, P, I, D, κ)
		tspan = (0.0, 100.0) 
		prob = ODEProblem(f, ψ, tspan) 
		sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)

		#Check for convergence 
		ψ = sol[:,end] 
		ω = ψ[1:n]
		θ = ψ[n+1:end]
		x = fsteadystate(θ, A, P, κ) 
		if norm(x,2) > 1e-3
			println("Warning: no steady state found.") 
		end

		#Compute the flow and normalise 
		flow_swing = edgepower(θ, E, κ) 
		flowmax_swing = maximum(abs.(flow_swing))
		fmaxes[z] = flowmax_swing
		flow_swing = flow_swing./flowmax_swing 

		redges = zeros(Float64, alphares)
		rtimes = zeros(Float64, alphares)

		for (ind, alpha) ∈ enumerate(alphas)

			#Reactive all the edges
			for i=1:length(net.edgelist)
				net.edgelist[i].active = true 
			end

			#Knock out most loaded edge 
			d = argmax(abs.(flow_swing))
			edgeset1 = filter(x->x≠d, edgeset) 
			E1 = E[1:end, 1:end .≠ d]
			net.edgelist[d].active = false 

			#Check to see if its broken 
			Adj = adjacencymatrix(E1)
			components, table = connectedcomponents(Adj)

			totedgessurviving = 0
			timetakenvec = []
			for i=1:components
				chunk = table[i] 
				n_subset = length(chunk)
				nodesubset = zeros(Int64, length(chunk))
				ω_subset = zeros(Float64, length(chunk))
				θ_subset = zeros(Float64, length(chunk))
				ψ_subset = zeros(Float64, 2*n_subset)
				for j=1:length(chunk)
					nodesubset[j] = nodeset[chunk[j]]
					ω_subset[j] = ω[chunk[j]]
					θ_subset[j] = θ[chunk[j]]
				end
				ψ_subset[1:n_subset] = ω_subset
				ψ_subset[n_subset+1:end] = θ_subset
				edgessurviving, timetaken = swingfracturewithtime!(net, ψ_subset, nodesubset, P, 3.0, alpha, κ, flowmax_swing, D)
				totedgessurviving += edgessurviving
				push!(timetakenvec, timetaken)
			end
			redges[ind] = totedgessurviving/m 
			rtimes[ind] = maximum(timetakenvec) 
		end

		#Find the critical value 
		aftertransition = findall(x->x>0.5, redges)
		if length(aftertransition) < 1 
			critvalues[z] = NaN 
		else
			ac = alphas[aftertransition[1]]
			critvalues[z] = ac 	
		end 

		#Store S 
		edgemat[z,:] = redges
		tmat[z,:] = rtimes
	end
	return alphas, mean(edgemat, dims=1)[1,:], mean(tmat, dims=1)[1,:] ,mean(fmaxes), critvalues  
end



function svt(ensemblesize::Int64, n::Int64, ns::Int64, nd::Int64, q::Float64, 
								k::Int64, alpha::Float64,
								I::Float64, D::Float64, κ::Float64)


	println("Starting")
	vedges = zeros(Float64, ensemblesize) 
	vtimes = zeros(Float64, ensemblesize)
	for z=1:ensemblesize
		println("Ensemble = ", z)

		net, E = smallworldnetwork(n, k, q)
		sources, sinks = sourcesinklocs(ns, nd, n)
		P = sourcesinkvector(sources, sinks, n, 1.0)
		m = numedges(net)
		nodeset = collect(1:n)
		edgeset = collect(1:m) 
		A = adjacencymatrix(E)
		ψ = rand(2*n) 

		#Time step to find steady state 
		f(u, p, t) = fswing(u, A, P, I, D, κ)
		tspan = (0.0, 100.0) 
		prob = ODEProblem(f, ψ, tspan) 
		sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)

		#Check for convergence 
		ψ = sol[:,end] 
		ω = ψ[1:n]
		θ = ψ[n+1:end]
		x = fsteadystate(θ, A, P, κ) 
		if norm(x,2) > 1e-3
			println("Warning: no steady state found.") 
		end

		#Compute the flow and normalise 
		flow_swing = edgepower(θ, E, κ) 
		flowmax_swing = maximum(abs.(flow_swing))
		flow_swing = flow_swing./flowmax_swing 


		#Reactive all the edges
		for i=1:length(net.edgelist)
			net.edgelist[i].active = true 
		end

		#Knock out most loaded edge 
		d = argmax(abs.(flow_swing))
		edgeset1 = filter(x->x≠d, edgeset) 
		E1 = E[1:end, 1:end .≠ d]
		net.edgelist[d].active = false 

		#Check to see if its broken 
		Adj = adjacencymatrix(E1)
		components, table = connectedcomponents(Adj)

		totedgessurviving = 0
		timetakenvec = []
		for i=1:components
			chunk = table[i] 
			n_subset = length(chunk)
			nodesubset = zeros(Int64, length(chunk))
			ω_subset = zeros(Float64, length(chunk))
			θ_subset = zeros(Float64, length(chunk))
			ψ_subset = zeros(Float64, 2*n_subset)
			for j=1:length(chunk)
				nodesubset[j] = nodeset[chunk[j]]
				ω_subset[j] = ω[chunk[j]]
				θ_subset[j] = θ[chunk[j]]
			end
			ψ_subset[1:n_subset] = ω_subset
			ψ_subset[n_subset+1:end] = θ_subset
			edgessurviving, timetaken = swingfracturewithtime!(net, ψ_subset, nodesubset, P, 3.0, alpha, κ, flowmax_swing, D)
			totedgessurviving += edgessurviving
			push!(timetakenvec, timetaken)
		end
		vedges[z] = totedgessurviving/m 
		vtimes[z] = maximum(timetakenvec) 
	end
	return vedges, vtimes   
end



function swingcascadebisection(ensemblesize::Int64, tol::Float64, n::Int64, ns::Int64, nd::Int64, 
								q::Float64, k::Int64, κ::Float64, I::Float64, D::Float64)


	#Make an ensemble of networks, with node placements and initial flow pattern
	ensemble = Array{CascadeNetworkSwing, 1}(undef, ensemblesize) 
	fmaxes = Array{Float64, 1}(undef, ensemblesize)
	for i=1:ensemblesize
		println(i)

		# Construct each network 
		net, E = smallworldnetwork(n, k, q)
		# net, E = austriannetwork() 
		sources, sinks = sourcesinklocsmod(net, ns, nd, n) 
		P = sourcesinkvector(sources, sinks, n, 1.0) 
		m = numedges(net) 

		#Make the initial state 
		A = adjacencymatrix(E) 
		ψ = rand(2*n) 
		f(u, p, t) = fswing(u, A, P, I, D, κ) 
		tspan = (0.0, 250.0) 
		prob = ODEProblem(f, ψ, tspan) 
		sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)

		#Check for convergence 
		ψ = sol[:,end] 
		ω = ψ[1:n] 
		θ = ψ[n+1:end] 
		x = fsteadystate(θ, A, P, κ) 
		if norm(x,2) > 1e-3
			println("Warning: no steady state found.")
		end 

		#Compute the flow pattern and normalise 
		flow_swing = edgepower(θ, E, κ) 
		flowmax_swing = maximum(abs.(flow_swing))
		fmaxes[i] = flowmax_swing 
		flow_swing = flow_swing./flowmax_swing 

		#Store the ensemble member 
		ensemblemember = CascadeNetworkSwing(net, E, sources, sinks, flow_swing, ψ) 
		ensemble[i] = ensemblemember 
	end

	alpha = 0.01
	stepsize = 0.3
	beneath0 = true 
	beneath1 = true 
	iterations = 1

	while abs(stepsize) > tol 

		# println(alpha) 
		runningsurvivors = Array{Float64, 1}(undef, ensemblesize) 
		for z=1:ensemblesize 

			#Unpack the network 
			net = ensemble[z].network
			E = ensemble[z].E 
			m = numedges(net) 
			sources = ensemble[z].sources 
			sinks = ensemble[z].sinks 
			F = ensemble[z].initialflow 
			ψ = ensemble[z].initialnodestate
			ω = ψ[1:n] 
			θ = ψ[n+1:end] 
			nodeset = collect(1:n) 
			edgeset = collect(1:m) 
			P = sourcesinkvector(sources, sinks, n, 1.0) 

			#Reactivate all edges 
			for j=1:length(net.edgelist)
				net.edgelist[j].active = true 
			end

			#Knock out largest
			d = argmax(abs.(F)) 
			edgeset1 = filter(x->x≠d, edgeset) 
			E1 = E[1:end, 1:end .≠ d]
			net.edgelist[d].active = false 

			#Check to see if its broken 
			Adj = adjacencymatrix(E1)
			components, table = connectedcomponents(Adj)

			#Do cascade
			totedgessurviving = 0
			for i=1:components
				chunk = table[i] 
				n_subset = length(chunk)
				nodesubset = zeros(Int64, length(chunk))
				ω_subset = zeros(Float64, length(chunk))
				θ_subset = zeros(Float64, length(chunk))
				ψ_subset = zeros(Float64, 2*n_subset)
				for j=1:length(chunk)
					nodesubset[j] = nodeset[chunk[j]]
					ω_subset[j] = ω[chunk[j]]
					θ_subset[j] = θ[chunk[j]]
				end
				ψ_subset[1:n_subset] = ω_subset
				ψ_subset[n_subset+1:end] = θ_subset
				edgessurviving = swingfracture!(net, ψ_subset, nodesubset, P, 3.0, alpha, κ, fmaxes[z])
				totedgessurviving += edgessurviving
			end

			#Record fraction of edges remaining
			runningsurvivors[z] = totedgessurviving/Float64(m)  
		end

		#If fraction has crossed the target then reverse and half the step size 
		global avremaining = mean(runningsurvivors) 
		# println("avremaining", avremaining)
		if avremaining > 0.5
			beneath1 = false 
		else 
			beneath1 = true
		end 
		if beneath1 ≠ beneath0 
			stepsize = -stepsize/2.0 
		end 
		beneath0 = beneath1 
		alpha += stepsize 
		iterations += 1 

	end
	return alpha, mean(fmaxes), iterations  
end


function criticalcouplingbisection(ensemblesize::Int64, tol::Float64, n::Int64, ns::Int64, nd::Int64, 
								q::Float64, k::Int64, κ::Float64, I::Float64, D::Float64)

	#Iterate over ensemble 
	crit = zeros(Float64, ensemblesize)
	κ0 = κ 
	for z=1:ensemblesize
		println(z)
		#Make the network 
		net, E = austriannetwork() 
		sources, sinks = sourcesinklocsmod(net, ns, nd, n) 
		P = sourcesinkvector(sources, sinks, n, 1.0) 
		m = numedges(net) 

		#Make the initial state 
		A = adjacencymatrix(E) 
		ψlast = rand(2*n) 
		ψnew = rand(2*n) 
		f(u, p, t) = fswing(u, A, P, I, D, κ) 
		tspan = (0.0, 200.0) 
		prob = ODEProblem(f, ψlast, tspan) 
		sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,
			dense=false,save_end=true) 
		ψlast = sol[:,end] 

		obj1 = NetStore2(A, P, E, κ, 0.0, 0.0, 0.0, true, false)

		tol = 1e-3
	
		#Iterate over coupling
		stepsize = 0.01 
		bisecting = true
		κ = κ0 
		κnew = κ
		κold = κ 
		while bisecting 
			println(κ)

			f(u, p, t) = fswing(u, A, P, I, D, κ) 
			obj1.kappa = κ 
			κnew = κ 
			tspan = (0.0, 100.0) 
			prob = ODEProblem(f, ψlast, tspan) 

			#Set up callback function 
			convergedflag = false
			function conditiontemp(u1, t1, integrator, ob::NetStore2)   
				θ = u1[n+1:end]
				if (t1 > 10.0)
					resid = fsteadystate(θ, ob.adjmat, ob.pvec, ob.kappa)
					if norm(resid,2) < 1e-5
						convergedflag = true 
						return true 
					end 
				end 
				return false 
			end
			condition(u,t,integrator) = conditiontemp(u,t,integrator,obj1) 
			affect!(integrator) = terminate!(integrator)
			cb = DiscreteCallback(condition, affect!) 
			sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,
				dense=false,save_end=true, callback=cb)  
			if convergedflag == true
				ψlast = sol[:,end]

				#If bisection done, then stop
				if stepsize < tol 
					crit[z] = κ 
					bisecting = false 
				end

				#Else, continue 
				κold = κ 
				κ = κ - stepsize 
			else 
				stepsize = stepsize/2.0
				κ = κold - stepsize 
			end
		end 
	end
	return crit

end


function swingcascadeouterbisection(ensemblesize::Int64, tol::Float64, n::Int64, ns::Int64, nd::Int64, 
								q::Float64, k::Int64, κ::Float64, I::Float64, D::Float64)

	
	rhos = zeros(Float64, ensemblesize)
	for z=1:ensemblesize 
		println("Ensemble ", z)
		net, E = smallworldnetwork(n, k, q) 
		sources, sinks = sourcesinklocsmod(net, ns, nd, n) 
		P = sourcesinkvector(sources, sinks, n, 1.0) 
		m = numedges(net)  

		#Make the initial state 
		A = adjacencymatrix(E) 
		ψ = rand(2*n) 
		f(u, p, t) = fswing(u, A, P, I, D, κ) 
		tspan = (0.0, 250.0) 
		prob = ODEProblem(f, ψ, tspan) 
		sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true) 

		#Check for convergence 
		alpha = 0.01
		ψ = sol[:,end] 
		ω = ψ[1:n] 
		θ = ψ[n+1:end] 
		x = fsteadystate(θ, A, P, κ) 
		if norm(x,2) > 1e-3
			println("Warning: no steady state found.")
			alpha = NaN 
		else  

			#Compute the flow pattern and normalise 
			F = edgepower(θ, E, κ) 
			fmax = maximum(abs.(F)) 
			F = F./fmax 

			#Now do bisection 
			stepsize = 0.3
			beneath0 = true 
			beneath1 = true 
			iterations = 1
			while abs(stepsize) > tol 

				#Reactivate all edges
				nodeset = collect(1:n) 
				edgeset = collect(1:m) 
				for j=1:length(net.edgelist)
					net.edgelist[j].active = true 
				end

				#Knock out largest
				d = argmax(abs.(F)) 
				edgeset1 = filter(x->x≠d, edgeset) 
				E1 = E[1:end, 1:end .≠ d]
				net.edgelist[d].active = false

				#Check to see if its broken 
				Adj = adjacencymatrix(E1)
				components, table = connectedcomponents(Adj)

				#Do cascade
				totedgessurviving = 0
				for i=1:components
					chunk = table[i] 
					n_subset = length(chunk)
					nodesubset = zeros(Int64, length(chunk))
					ω_subset = zeros(Float64, length(chunk))
					θ_subset = zeros(Float64, length(chunk))
					ψ_subset = zeros(Float64, 2*n_subset)
					for j=1:length(chunk)
						nodesubset[j] = nodeset[chunk[j]]
						ω_subset[j] = ω[chunk[j]]
						θ_subset[j] = θ[chunk[j]]
					end
					ψ_subset[1:n_subset] = ω_subset
					ψ_subset[n_subset+1:end] = θ_subset
					edgessurviving = swingfracture!(net, ψ_subset, nodesubset, P, 2.0, alpha, κ, fmax)
					totedgessurviving += edgessurviving
				end
				S = totedgessurviving/Float64(m)

				#If fraction has crossed the target then reverse and half the step size 
				if S > 0.5
					beneath1 = false 
				else 
					beneath1 = true
				end 
				if beneath1 ≠ beneath0 
					stepsize = -stepsize/2.0 
				end 
				beneath0 = beneath1 
				alpha += stepsize 
				iterations += 1 

			end 
		end 
		rhos[z] = alpha 
	end 
	return rhos 

end


function swingcascadebreakdownvsalpha(ensemblesize::Int64, n::Int64, ns::Int64, nd::Int64, q::Float64, 
								k::Int64, alphares::Int64, alphamin::Float64, alphamax::Float64,
								I::Float64, D::Float64, κ::Float64)

	alphas = collect(range(alphamin, length=alphares, stop=alphamax))
	edgemat = zeros(Float64, ensemblesize, alphares) 
	overloadmat = zeros(Float64, ensemblesize, alphares)
	desyncmat = zeros(Float64, ensemblesize, alphares)
	unbalancemat = zeros(Float64, ensemblesize, alphares)
	fmaxes = zeros(Float64, ensemblesize)
	critvalues = zeros(Float64, ensemblesize) 

	println("Starting")
	for z=1:ensemblesize
		println("Ensemble = ", z)

		net, E = smallworldnetwork(n, k, q)
		# net, E = scalefreenetwork(n, 2*n, 3) 
		sources, sinks = sourcesinklocs(ns, nd, n)
		P = sourcesinkvector(sources, sinks, n, 1.0)
		m = numedges(net)
		nodeset = collect(1:n)
		edgeset = collect(1:m) 
		A = adjacencymatrix(E)
		ψ = rand(2*n) 

		#Time step to find steady state 
		f(u, p, t) = fswing(u, A, P, I, D, κ)
		tspan = (0.0, 100.0) 
		prob = ODEProblem(f, ψ, tspan) 
		sol = solve(prob,reltol=1e-8,abstol=1e-8,save_everystep=false,dense=false,save_end=true)

		#Check for convergence 
		ψ = sol[:,end] 
		ω = ψ[1:n]
		θ = ψ[n+1:end]
		x = fsteadystate(θ, A, P, κ) 
		if norm(x,2) > 1e-3
			println("Warning: no steady state found.") 
		end

		#Compute the flow and normalise 
		flow_swing = edgepower(θ, E, κ) 
		flowmax_swing = maximum(abs.(flow_swing))
		fmaxes[z] = flowmax_swing
		flow_swing = flow_swing./flowmax_swing 

		redges = zeros(Float64, alphares)
		roverloads = zeros(Float64, alphares)
		rdesyncs = zeros(Float64, alphares)
		runbalances = zeros(Float64, alphares)

		for (ind, alpha) ∈ enumerate(alphas)

			#Reactive all the edges
			for i=1:length(net.edgelist)
				net.edgelist[i].active = true 
			end

			#Knock out most loaded edge 
			d = argmax(abs.(flow_swing))
			edgeset1 = filter(x->x≠d, edgeset) 
			E1 = E[1:end, 1:end .≠ d]
			net.edgelist[d].active = false 

			#Check to see if its broken 
			Adj = adjacencymatrix(E1)
			components, table = connectedcomponents(Adj)

			totedgessurviving = 0
			totaloverloads = 0 
			totaldesyncs = 0 
			totalunbalances = 0 
			for i=1:components
				chunk = table[i] 
				n_subset = length(chunk)
				nodesubset = zeros(Int64, length(chunk))
				ω_subset = zeros(Float64, length(chunk))
				θ_subset = zeros(Float64, length(chunk))
				ψ_subset = zeros(Float64, 2*n_subset)
				for j=1:length(chunk)
					nodesubset[j] = nodeset[chunk[j]]
					ω_subset[j] = ω[chunk[j]]
					θ_subset[j] = θ[chunk[j]]
				end
				ψ_subset[1:n_subset] = ω_subset
				ψ_subset[n_subset+1:end] = θ_subset
				edgessurviving, overloads, desyncs, unbalances = swingfracturewithbreakdown!(net, ψ_subset, nodesubset, P, 3.0, alpha, κ, flowmax_swing)
				totedgessurviving += edgessurviving
				totaloverloads += overloads 
				totaldesyncs += desyncs 
				totalunbalances += unbalances 
			end
			redges[ind] = totedgessurviving/m 
			roverloads[ind] = (m - totedgessurviving - totaldesyncs - totalunbalances)/m 
			rdesyncs[ind] = totaldesyncs/m
			runbalances[ind] = totalunbalances/m
		end

		#Find the critical value 
		aftertransition = findall(x->x>0.5, redges)
		if length(aftertransition) < 1 
			critvalues[z] = NaN 
		else
			ac = alphas[aftertransition[1]]
			critvalues[z] = ac 	
		end 

		#Store S 
		edgemat[z,:] = redges
		overloadmat[z,:] = roverloads
		desyncmat[z,:] = rdesyncs
		unbalancemat[z,:] = runbalances
	end
	return alphas, mean(edgemat, dims=1)[1,:], mean(overloadmat, dims=1)[1,:], mean(desyncmat, dims=1)[1,:], mean(unbalancemat, dims=1)[1,:], mean(fmaxes), critvalues  
end






















