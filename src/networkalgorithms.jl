function depthfirstsearch!(A::Matrix{Float64}, u::Int64, 
			visited::Vector{Bool}, contents::Vector{Int64})
	n = size(A,1)
	for v=1:n
		adj = A[u, v] 
		if (abs(adj) > 0.5) && (u != v)
			if visited[v] == false
				visited[v] = true 
				push!(contents, v)
				contents = depthfirstsearch!(A, v, visited, contents)
			end 
		end 
	end
	return contents 
end


function connectedcomponents(A::Matrix{Float64})
	n = size(A,1) 
	components = 0 
	visited = Array{Bool, 1}(undef, n)
	visited[:] .= false 
	table = Array{Array{Int64, 1}, 1}([]) 
	for u=1:n
		if visited[u] == false
			visited[u] = true
			components += 1 
			contents = Array{Int64, 1}([]) 
			push!(contents, u)
			contents = depthfirstsearch!(A, u, visited, contents) 
			push!(table, contents)
		end
	end
	return components, table  
end


function partitionadjacency(A::Matrix{Float64}, v::Vector{Int64})
	n = size(A,1)
	n1 = size(v,1)
	Atemp = zeros(Float64, n, n1)
	Afinal = zeros(Float64, n1, n1) 
	for (i, val) ∈ enumerate(v)
		Atemp[:,i] = A[:,val] 
	end
	for (i, val) ∈ enumerate(v)
		Afinal[i,:] = Atemp[val,:]
	end
	return Afinal 
end


function avcompontentsize(A::Matrix{Float64})
	components, table = connectedcomponents(A) 
	s = 0
	for i=1:components
		s += size(table[i],1)
	end
	return s/components 
end


function clustersize(A::Matrix{Float64}, v::Vector{Int64})
	Anew = partitionadjacency(A, v) 
	avs = avcompontentsize(Anew)
	return avs 
end



function nodevoltages(L::Matrix{Float64}, P::Vector{Float64})
	n = size(L,1)
	Lred = L[2:end,2:end]
	C = zeros(Float64, n, n)
	C[2:end,2:end] = inv(Lred)
	return C*P 
end 


function edgecurrents(theta::Vector{Float64}, E::Matrix{Float64}, X::Vector{Float64})
	tmp = transpose(E)*theta
	return tmp./X
end


function fracture!(net::Graph, nodeset::Vector{Int64}, P::Vector{Float64}, X::Vector{Float64}, Finit::Vector{Float64}, info::NetworkRecorder, t::Int64, alpha::Float64, fmaxtemp::Float64, seq::Bool)

	tol = 1e-5 
	simtime = t + 1

	#Check the status of P, see if has any sources or sinks. If not, halt recursion
	nodecardinality = length(nodeset) 
	P1 = zeros(Float64, nodecardinality) 
	sourcecounter = 0
	weighedsourcecounter = 0.0 
	sinkcounter = 0 
	for i=1:nodecardinality
		vi = nodeset[i] 
		P1[i] = P[vi] 
		if P1[i] > tol
			sourcecounter += 1
			weighedsourcecounter += P1[i]
		end
		if P1[i] < -tol
			sinkcounter += 1
		end 
	end 
	if (sourcecounter == 0) || (sinkcounter == 0)
		if simtime ∉ info.simulationtime
			push!(info.simulationtime, simtime) 
			push!(info.liveedges, 0)
			push!(info.livesources, 0)
			push!(info.weightedlivesources, 0.0) 
			push!(info.livesinks, 0)
			push!(info.totalpower, 0.0)
			if seq
				flowvec = zeros(Float64, info.liveedges[1]) 
				flowvec[:] .= NaN 
				push!(info.flowsequence, flowvec)
			end 
		end
		return info 
	end

	# #Balance P1 - if there is power excess, increase loads
	# delta = sum(P1) 
	# individualdelta = delta/Float64(sinkcounter)
	# for i=1:nodecardinality
	# 	if P1[i] < -tol 
	# 		P1[i] -= individualdelta
	# 	end 
	# end   

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

	# # Balance P1 heterogeneously  
	# delta = sum(P1)/2.0 #This is the amount of surplus given to sources (or sinks) to deal with
	# avP = 0.0 
	# for i=1:nodecardinality
	# 	if P1[i] > tol 
	# 		avP += P1[i]
	# 	end
	# end 
	# avP /= sourcecounter
	# for i=1:nodecardinality
	# 	if P1[i] < -tol 
	# 		P1[i] -=  abs(P1[i])*delta/(Float64(sinkcounter)*avP)  
	# 	end 
	# 	if P1[i] > tol
	# 		P1[i] -= abs(P1[i])*delta/(Float64(sourcecounter)*avP)
	# 	end
	# end  
	# println("sum", "    ", sum(P1))


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

	#Construct X
	X1 = zeros(Float64, edgecardinality)
	for i=1:edgecardinality
		X1[i] = X[edgeset[i]]
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

	# A = adjacencymatrix(E1)
	# comt, tab = connectedcomponents(A)
	# println("com... = ", comt)
	#Solve the DC flow equations 
	L = weightedlaplacianmatrix(E1, 1.0./X1) 
	theta = nodevoltages(L, P1)
	F = edgecurrents(theta, E1, X1)./fmaxtemp

	#Record the number of edges in a tuple of recursion depth and number of edges 
	if simtime ∈ info.simulationtime
		i = findall(y->y==simtime, info.simulationtime) 
		info.liveedges[i] .+= edgecardinality 
		info.livesources[i] .+= sourcecounter
		info.weightedlivesources .+= weighedsourcecounter
		info.livesinks[i] .+= sinkcounter
		info.totalpower[i] .+= sum(abs.(P1))/2.0
		if seq 
			for j=1:edgecardinality
				ej = edgeset[j]
				info.flowsequence[i][ej] = abs(F[j])  
			end
		end 
	else
		push!(info.simulationtime, simtime)
		push!(info.liveedges, edgecardinality)
		push!(info.livesources, sourcecounter)
		push!(info.weightedlivesources, weighedsourcecounter) 
		push!(info.livesinks, sinkcounter)
		push!(info.totalpower, sum(abs.(P1))/2.0)
		if seq
			flowvec = zeros(Float64, info.liveedges[1]) 
			flowvec[:] .= NaN 
			for j=1:edgecardinality
				ej = edgeset[j]
				flowvec[ej] = abs(F[j])
			end
			push!(info.flowsequence, flowvec)
		end
	end

	#Scan through new flow vector and record any overloads
	#If no new overloads, halt recursion
	overloadcounter = 0
	survivors = Array{Int64, 1}([]) 
	for i=1:edgecardinality
		ei = edgeset[i]
		if abs(F[i]) > alpha #*abs(Finit[ei]) 
			net.edgelist[ei].active = false
			overloadcounter += 1 
		else
			push!(survivors, i)
		end 
	end
	if overloadcounter == 0 
		return info 
	end 

	#Delete overloads from E1 
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
	L = nothing
	F = nothing

	#Call fracture function on each of the components, passing down info
	for i=1:components
		chunk = table[i] 
		nodesubset = zeros(Int64, length(chunk))
		for j=1:length(chunk)
			nodesubset[j] = nodeset[chunk[j]] 
		end
		info = fracture!(net, nodesubset, P, X, Finit, info, simtime, alpha, fmaxtemp, seq) 
	end 

	return info 
end



function fracture2!(net::Graph, nodeset::Vector{Int64}, P::Vector{Float64}, X::Vector{Float64}, Finit::Vector{Float64}, t::Int64, alpha::Float64, fmaxtemp::Float64)

	tol = 1e-5 
	simtime = t + 1

	#Check the status of P, see if has any sources or sinks. If not, halt recursion
	nodecardinality = length(nodeset) 
	P1 = zeros(Float64, nodecardinality) 
	sourcecounter = 0
	weighedsourcecounter = 0.0 
	sinkcounter = 0 
	for i=1:nodecardinality
		vi = nodeset[i] 
		P1[i] = P[vi] 
		if P1[i] > tol
			sourcecounter += 1
			weighedsourcecounter += P1[i]
		end
		if P1[i] < -tol
			sinkcounter += 1
		end 
	end 
	if (sourcecounter == 0) || (sinkcounter == 0)
		return 0 
	end

	# #Balance P1 - if there is power excess, increase loads
	# delta = sum(P1) 
	# individualdelta = delta/Float64(sinkcounter)
	# for i=1:nodecardinality
	# 	if P1[i] < -tol 
	# 		P1[i] -= individualdelta
	# 	end 
	# end  

	# #Balance P1 homogeneously 
	# delta = sum(P1)/2.0 
	# individualdelta_sink = delta/Float64(sinkcounter)
	# individualdelta_source = delta/Float64(sourcecounter)
	# for i=1:nodecardinality
	# 	if P1[i] < -tol 
	# 		P1[i] -= individualdelta_sink  
	# 	end 
	# 	if P1[i] > tol
	# 		P1[i] -= individualdelta_source
	# 	end
	# end 

	# Balance P1 heterogeneously  
	delta = sum(P1)/2.0 #This is the amount of surplus given to sources (or sinks) to deal with
	avP = 0.0 
	for i=1:nodecardinality
		if P1[i] > tol 
			avP += P1[i]
		end
	end 
	avP /= sourcecounter
	for i=1:nodecardinality
		if P1[i] < -tol 
			P1[i] -=  abs(P1[i])*delta/(Float64(sinkcounter)*avP)  
		end 
		if P1[i] > tol
			P1[i] -= abs(P1[i])*delta/(Float64(sourcecounter)*avP)
		end
	end  
	# println("sum", "    ", sum(P1))


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

	#Construct X
	X1 = zeros(Float64, edgecardinality)
	for i=1:edgecardinality
		X1[i] = X[edgeset[i]]
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

	#Solve the DC flow equations 
	L = weightedlaplacianmatrix(E1, 1.0./X1) 
	theta = nodevoltages(L, P1)
	F = edgecurrents(theta, E1, X1)./fmaxtemp

	#Scan through new flow vector and record any overloads
	overloadcounter = 0
	survivors = Array{Int64, 1}([]) 
	for i=1:edgecardinality
		ei = edgeset[i]
		if abs(F[i]) > alpha #*abs(Finit[ei]) 
			net.edgelist[ei].active = false
			overloadcounter += 1 
		else
			push!(survivors, i)
		end 
	end

	#If no new overloads, halt recursion
	if overloadcounter == 0 
		return edgecardinality  
	end 

	#Delete overloads from E1 
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
	L = nothing
	F = nothing

	#Call fracture function on each of the components, passing down info
	descendent_survivingedges = 0 
	for i=1:components
		chunk = table[i] 
		nodesubset = zeros(Int64, length(chunk))
		for j=1:length(chunk)
			nodesubset[j] = nodeset[chunk[j]] 
		end
		survivingedges = fracture2!(net, nodesubset, P, X, Finit, simtime, alpha, fmaxtemp) 
		descendent_survivingedges += survivingedges
	end 

	return descendent_survivingedges  
end


function fracturerand!(net::Graph, nodeset::Vector{Int64}, P::Vector{Float64}, X::Vector{Float64}, Finit::Vector{Float64}, capnoise::Vector{Float64}, t::Int64, alpha::Float64, fmaxtemp::Float64)

	tol = 1e-5 
	simtime = t + 1

	#Check the status of P, see if has any sources or sinks. If not, halt recursion
	nodecardinality = length(nodeset) 
	P1 = zeros(Float64, nodecardinality) 
	sourcecounter = 0
	weighedsourcecounter = 0.0 
	sinkcounter = 0 
	for i=1:nodecardinality
		vi = nodeset[i] 
		P1[i] = P[vi] 
		if P1[i] > tol
			sourcecounter += 1
			weighedsourcecounter += P1[i]
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

	#Construct X
	X1 = zeros(Float64, edgecardinality)
	for i=1:edgecardinality
		X1[i] = X[edgeset[i]]
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

	#Solve the DC flow equations 
	L = weightedlaplacianmatrix(E1, 1.0./X1) 
	theta = nodevoltages(L, P1)
	F = edgecurrents(theta, E1, X1)./fmaxtemp

	#Scan through new flow vector and record any overloads
	overloadcounter = 0
	survivors = Array{Int64, 1}([]) 
	for i=1:edgecardinality
		ei = edgeset[i]
		if abs(F[i]) > alpha*(1.0 + capnoise[ei])    
			net.edgelist[ei].active = false
			overloadcounter += 1 
		else
			push!(survivors, i)
		end 
	end

	#If no new overloads, halt recursion
	if overloadcounter == 0 
		return edgecardinality  
	end 

	#Delete overloads from E1 
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
	L = nothing
	F = nothing

	#Call fracture function on each of the components, passing down info
	descendent_survivingedges = 0 
	for i=1:components
		chunk = table[i] 
		nodesubset = zeros(Int64, length(chunk))
		for j=1:length(chunk)
			nodesubset[j] = nodeset[chunk[j]] 
		end
		survivingedges = fracture2!(net, nodesubset, P, X, Finit, simtime, alpha, fmaxtemp) 
		descendent_survivingedges += survivingedges
	end 

	return descendent_survivingedges  
end






function fracture3!(net::Graph, nodeset::Vector{Int64}, P::Vector{Float64}, X::Vector{Float64}, Finit::Vector{Float64}, t::Int64, alpha::Float64, fmaxtemp::Float64)

	tol = 1e-5 
	simtime = t + 1

	#Check the status of P, see if has any sources or sinks. If not, halt recursion
	nodecardinality = length(nodeset) 
	P1 = zeros(Float64, nodecardinality) 
	sourcecounter = 0
	weighedsourcecounter = 0.0 
	sinkcounter = 0 
	for i=1:nodecardinality
		vi = nodeset[i] 
		P1[i] = P[vi] 
		if P1[i] > tol
			sourcecounter += 1
			weighedsourcecounter += P1[i]
		end
		if P1[i] < -tol
			sinkcounter += 1
		end 
	end 
	if (sourcecounter == 0) || (sinkcounter == 0)
		return 1 
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

	#Construct X
	X1 = zeros(Float64, edgecardinality)
	for i=1:edgecardinality
		X1[i] = X[edgeset[i]]
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

	#Solve the DC flow equations 
	L = weightedlaplacianmatrix(E1, 1.0./X1) 
	theta = nodevoltages(L, P1)
	F = edgecurrents(theta, E1, X1)./fmaxtemp

	#Scan through new flow vector and record any overloads
	overloadcounter = 0
	survivors = Array{Int64, 1}([]) 
	for i=1:edgecardinality
		ei = edgeset[i]
		if abs(F[i]) > alpha #*abs(Finit[ei]) 
			net.edgelist[ei].active = false
			overloadcounter += 1 
		else
			push!(survivors, i)
		end 
	end

	#If no new overloads, halt recursion
	if overloadcounter == 0 
		return 1  
	end 

	#Delete overloads from E1 
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
	L = nothing
	F = nothing

	#Call fracture function on each of the components, passing down info
	descendent_leaves = 0 
	for i=1:components
		chunk = table[i] 
		nodesubset = zeros(Int64, length(chunk))
		for j=1:length(chunk)
			nodesubset[j] = nodeset[chunk[j]] 
		end
		leaves = fracture2!(net, nodesubset, P, X, Finit, simtime, alpha, fmaxtemp) 
		descendent_leaves += leaves
	end 

	return descendent_leaves 
end
















