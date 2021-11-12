using .NetworkResilience
using PyPlot

ac, fmax, iter = cascadebisection(3000, 1e-4, 100, 50, 50, 0.0, 4, 1.0)
println(iter)

# alpha = 10.0
# net, E0, X, P = testnetwork() 
# n = numnodes(net)
# m = numedges(net)
# nodeset = collect(1:n)
# edgeset0 = collect(1:m) 

# #Solve the DC flow equations
# L = weightedlaplacianmatrix(E0, 1.0./X)  
# theta = nodevoltages(L, P)
# F = edgecurrents(theta, E0, X)

# #Knock out an edge  
# edgeset1 = filter(x->x≠2, edgeset0) 
# E1 = E0[1:end, 1:end .≠ 2]
# net.edgelist[2].active = false 

# #Check to see if its broken 
# Adj = adjacencymatrix(E1)
# components, table = connectedcomponents(Adj) 

# #stuff needed for function call 
# info = [1 length(edgeset1)+1] 
# nodeset = table[1] 
# simtime = 1

# ###############################
# tol = 1e-5 
# simtime += 1

# #Check the status of P, see if has any sources or sinks. If not, halt recursion
# nodecardinality = length(nodeset)
# P1 = zeros(Float64, nodecardinality) 
# sourcecounter = 0
# sinkcounter = 0 
# for i=1:nodecardinality
# 	vi = nodeset[i] 
# 	P1[i] = P[vi]
# 	if P1[i] > tol
# 		global sourcecounter += 1
# 	end
# 	if P1[i] < -tol
# 		global sinkcounter += 1
# 	end 
# end 
# if (sourcecounter == 0) || (sinkcounter == 0)
# 	println("poop")
# end

# #Balance P1 - if there is power excess, increase loads
# delta = sum(P1) 
# individualdelta = delta/Float64(sinkcounter)
# for i=1:nodecardinality
# 	if P1[i] < -tol 
# 		P1[i] -= individualdelta
# 	end 
# end  

# #Find which edges are connected to the node-set
# edgeset = Array{Int64, 1}([]) 
# for vi ∈ nodeset
# 	v = net.nodelist[vi]
# 	for e ∈ v.outedges
# 		if (e ∉ edgeset) && (net.edgelist[e].active == true)
# 			push!(edgeset, e)
# 		end 	
# 	end
# 	for e ∈ v.inedges 
# 		if (e ∉ edgeset) && (net.edgelist[e].active == true)
# 			push!(edgeset, e)
# 		end 	
# 	end 
# end
# edgecardinality = length(edgeset) 
# # println(edgeset)

# #Construct X
# # println(length(edgeset)) 
# X1 = zeros(Float64, edgecardinality)
# for i=1:edgecardinality
# 	X1[i] = X[edgeset[i]]
# end  

# #Record the number of edges in a tuple of recursion depth and number of edges 
# if simtime ∈ info[:,1]
# 	i = findall(y->y==simtime, info[:,1])
# 	info[i,2] += edgecardinality
# else
# 	tmp = [simtime edgecardinality]
# 	info = vcat(info, tmp)
# end

# #Construct E1 
# E1 = zeros(Float64, nodecardinality, edgecardinality) 
# for i=1:nodecardinality
# 	vi = nodeset[i] 
# 	v = net.nodelist[vi]
# 	for j=1:edgecardinality
# 		if edgeset[j] ∈ v.outedges
# 			E1[i,j] = 1.0
# 		end
# 		if edgeset[j] ∈ v.inedges
# 			E1[i,j] = -1.0
# 		end 
# 	end		
# end

# #Solve the DC flow equations 
# L = weightedlaplacianmatrix(E1, 1.0./X1) 
# theta = nodevoltages(L, P1)
# F1 = edgecurrents(theta, E1, X1)

# # Scan through new flow vector and record any overloads
# # If no new overloads, halt recursion
# overloadcounter = 0
# survivors = Array{Int64, 1}([]) 
# println(edgeset)
# for i=1:edgecardinality
# 	ei = edgeset[i]
# 	if abs(F1[i]) > (1.0+alpha)*abs(F[ei]) 
# 		net.edgelist[ei].active = false
# 		global overloadcounter += 1 
# 	else
# 		push!(survivors, i)
# 		println(ei)
# 	end 
# end
# if overloadcounter == 0
# 	println("poop")
# end 
# # survivors = Array{Int64, 1}([])
# # for i=1:edgecardinality
# # 	ei = edgeset[i] 
# # 	if ei == 5
# # 		net.edgelist[ei].active = false
# # 	else
# # 		push!(survivors, i) 
# # 	end
# # end

# #Delete overloads from E1 
# newedgecard = length(survivors) 
# E2 = zeros(Float64, nodecardinality, newedgecard) 
# for j=1:newedgecard
# 	E2[:,j] = E1[:,survivors[j]] 
# end

# #Check to see if its broken 
# Adj = adjacencymatrix(E2)
# components, table = connectedcomponents(Adj)

# println(info)

# ############################################################
# nodeset = table[1]
# simtime += 1 

# #Check the status of P, see if has any sources or sinks. If not, halt recursion
# nodecardinality = length(nodeset)
# P1 = zeros(Float64, nodecardinality) 
# sourcecounter = 0
# sinkcounter = 0 
# for i=1:nodecardinality
# 	vi = nodeset[i] 
# 	P1[i] = P[vi]
# 	if P1[i] > tol
# 		global sourcecounter += 1
# 	end
# 	if P1[i] < -tol
# 		global sinkcounter += 1
# 	end 
# end 
# if (sourcecounter == 0) || (sinkcounter == 0)
# 	println("poop")
# end

# #Balance P1 - if there is power excess, increase loads
# delta = sum(P1) 
# individualdelta = delta/Float64(sinkcounter)
# for i=1:nodecardinality
# 	if P1[i] < -tol 
# 		P1[i] -= individualdelta
# 	end 
# end 
# println(P1)

# #Find which edges are connected to the node-set
# edgeset = Array{Int64, 1}([]) 
# for vi ∈ nodeset
# 	v = net.nodelist[vi]
# 	for e ∈ v.outedges
# 		if (e ∉ edgeset) && (net.edgelist[e].active == true)
# 			push!(edgeset, e)
# 		end 	
# 	end
# 	for e ∈ v.inedges 
# 		if (e ∉ edgeset) && (net.edgelist[e].active == true)
# 			push!(edgeset, e)
# 		end 	
# 	end 
# end
# edgecardinality = length(edgeset)

# #Construct X
# println(length(edgeset)) 
# X1 = zeros(Float64, edgecardinality)
# for i=1:edgecardinality
# 	X1[i] = X[edgeset[i]]
# end  

# #Record the number of edges in a tuple of recursion depth and number of edges 
# if simtime ∈ info[:,1]
# 	i = findall(y->y==simtime, info[:,1])
# 	info[i,2] += edgecardinality
# else
# 	tmp = [simtime edgecardinality]
# 	info = vcat(info, tmp)
# end

# #Construct E1 
# E1 = zeros(Float64, nodecardinality, edgecardinality) 
# for i=1:nodecardinality
# 	vi = nodeset[i] 
# 	v = net.nodelist[vi]
# 	for j=1:edgecardinality
# 		if edgeset[j] ∈ v.outedges
# 			E1[i,j] = 1.0
# 		end
# 		if edgeset[j] ∈ v.inedges
# 			E1[i,j] = -1.0
# 		end 
# 	end		
# end

# #Solve the DC flow equations 
# L = weightedlaplacianmatrix(E1, 1.0./X1) 
# theta = nodevoltages(L, P1)
# F = edgecurrents(theta, E1, X1)

# overloadcounter = 0
# survivors = Array{Int64, 1}([]) 
# println(edgeset)
# for i=1:edgecardinality
# 	ei = edgeset[i]
# 	if abs(F1[i]) > (1.0+alpha)*abs(F[ei]) 
# 		net.edgelist[ei].active = false
# 		global overloadcounter += 1 
# 	else
# 		push!(survivors, i)
# 		println(ei)
# 	end 
# end
# if overloadcounter == 0
# 	println("poop")
# end 

# #Delete overloads from E1 
# newedgecard = length(survivors) 
# E2 = zeros(Float64, nodecardinality, newedgecard) 
# for j=1:newedgecard
# 	E2[:,j] = E1[:,survivors[j]] 
# end

# #Check to see if its broken 
# Adj = adjacencymatrix(E2)
# components, table = connectedcomponents(Adj)

# println(info)

############################################################
# nodeset = table[1]
# simtime += 1 

# #Check the status of P, see if has any sources or sinks. If not, halt recursion
# nodecardinality = length(nodeset)
# P1 = zeros(Float64, nodecardinality) 
# sourcecounter = 0
# sinkcounter = 0 
# for i=1:nodecardinality
# 	vi = nodeset[i] 
# 	P1[i] = P[vi]
# 	if P1[i] > tol
# 		global sourcecounter += 1
# 	end
# 	if P1[i] < -tol
# 		global sinkcounter += 1
# 	end 
# end 
# if (sourcecounter == 0) || (sinkcounter == 0)
# 	println("ending")
# end

# #Balance P1 - if there is power excess, increase loads
# delta = sum(P1) 
# individualdelta = delta/Float64(sinkcounter)
# for i=1:nodecardinality
# 	if P1[i] < -tol 
# 		P1[i] -= individualdelta
# 	end 
# end 
# println(P1)

# #Find which edges are connected to the node-set
# edgeset = Array{Int64, 1}([]) 
# for vi ∈ nodeset
# 	v = net.nodelist[vi]
# 	for e ∈ v.outedges
# 		if (e ∉ edgeset) && (net.edgelist[e].active == true)
# 			push!(edgeset, e)
# 		end 	
# 	end
# 	for e ∈ v.inedges 
# 		if (e ∉ edgeset) && (net.edgelist[e].active == true)
# 			push!(edgeset, e)
# 		end 	
# 	end 
# end
# edgecardinality = length(edgeset)

# #Construct X
# X1 = zeros(Float64, edgecardinality)
# for i=1:edgecardinality
# 	X1[i] = X[edgeset[i]]
# end  

# #Record the number of edges in a tuple of recursion depth and number of edges 
# if simtime ∈ info[:,1]
# 	i = findall(y->y==simtime, info[:,1])
# 	info[i,2] += edgecardinality
# else
# 	tmp = [simtime edgecardinality]
# 	info = vcat(info, tmp)
# end

# #Construct E1 
# E1 = zeros(Float64, nodecardinality, edgecardinality) 
# for i=1:nodecardinality
# 	vi = nodeset[i] 
# 	v = net.nodelist[vi]
# 	for j=1:edgecardinality
# 		if edgeset[j] ∈ v.outedges
# 			E1[i,j] = 1.0
# 		end
# 		if edgeset[j] ∈ v.inedges
# 			E1[i,j] = -1.0
# 		end 
# 	end		
# end

# #Solve the DC flow equations 
# L = weightedlaplacianmatrix(E1, 1.0./X1) 
# theta = nodevoltages(L, P1)
# F = edgecurrents(theta, E1, X1)

# survivors = Array{Int64, 1}([])
# for i=1:edgecardinality
# 	ei = edgeset[i] 
# 	if ei == 4
# 		net.edgelist[ei].active = false
# 	else
# 		push!(survivors, i) 
# 	end
# end

# #Delete overloads from E1 
# newedgecard = length(survivors) 
# E2 = zeros(Float64, nodecardinality, newedgecard) 
# for j=1:newedgecard
# 	E2[:,j] = E1[:,survivors[j]] 
# end

# #Check to see if its broken 
# Adj = adjacencymatrix(E2)
# components, table = connectedcomponents(Adj)











