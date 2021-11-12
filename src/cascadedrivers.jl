mutable struct NetworkRecorder
	simulationtime::Array{Int64,1}
	liveedges::Array{Int64,1}
	livesources::Array{Int64,1}
	weightedlivesources::Array{Float64,1}
	livesinks::Array{Int64,1} 
	totalpower::Array{Float64,1} 
	flowsequence::Array{Array{Float64,1},1}
end


function newrecorder(numedges::Int64, ns::Int64, P::Float64, nd::Int64, F::Float64, flowvec::Vector{Float64})
	r = NetworkRecorder(Int64[],Int64[],Int64[],Float64[],Int64[],Float64[],Float64[])
	push!(r.simulationtime, 0)
	push!(r.liveedges, numedges)
	push!(r.livesources, ns)
	push!(r.weightedlivesources, P) 
	push!(r.livesinks, nd) 
	push!(r.totalpower, F)
	push!(r.flowsequence, flowvec)
	return r 
end


mutable struct NetworkRecorder2
	iterations::Array{Int64,1}
	liveedges::Array{Int64,1}
	livesources::Array{Int64,1}
	livesinks::Array{Int64,1} 
	totalpower::Array{Float64,1} 
end


function cascaderesults(info::NetworkRecorder)
	timetaken = info.simulationtime[end]
	fracedges = Float64(info.liveedges[end])/Float64(info.liveedges[1])
	fracsources = Float64(info.livesources[end])/Float64(info.livesources[1])
	fracweightedsources = info.weightedlivesources[end]/info.weightedlivesources[1] 
	fracsinks = Float64(info.livesinks[end])/Float64(info.livesinks[1])
	fracpower = info.totalpower[end]/info.totalpower[1] 
	return timetaken, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, info.flowsequence
end


mutable struct CascadeNetwork
	network::Graph 
	E::Matrix{Float64}
	sources::Vector{Int64}
	sinks::Vector{Int64} 
	initialflow::Vector{Float64} 
end


mutable struct CascadeNetworkSwing
	network::Graph 
	E::Matrix{Float64}
	sources::Vector{Int64}
	sinks::Vector{Int64} 
	initialflow::Vector{Float64} 
	initialnodestate::Vector{Float64} 
end


mutable struct CascadeNetworkRand
	network::Graph 
	E::Matrix{Float64}
	sources::Vector{Int64}
	sinks::Vector{Int64} 
	initialflow::Vector{Float64}
	alphanoise::Vector{Float64} 
end


function triggerhist(F::Vector{Float64})
	m = length(F) 
	println(sum(abs.(F)))
	probs = Array{Float64, 1}(undef, m) 
	for i=1:m
		probs[i] = abs(F[i])/sum(abs.(F)) 
	end
	return probs 
end


function triggerindex(F::Vector{Float64})
	m = length(F)
	while true 
		z = rand(1:m)
		if rand() < abs(F[z])/sum(abs.(F))
			return z 
		end
	end 
end 


function cascadebisectionrandomtrigger(supensemblesize::Int64, subensemblesize::Int64, tol::Float64, n::Int64, ns::Int64, nd::Int64, q::Float64, k::Int64, Ptot::Float64) 


	#Make an ensemble of networks here, with node placements
	#and with initial flow pattern
	ensemble = Array{CascadeNetwork, 1}(undef, supensemblesize*subensemblesize)
	fmaxes = Array{Float64, 1}(undef, supensemblesize*subensemblesize) 
	degmeans = Array{Float64, 1}(undef, supensemblesize*subensemblesize) 
	triggervec = Array{Int64, 1}(undef, supensemblesize*subensemblesize) 
	for i=1:supensemblesize

		#Construct each network
		net, E = smallworldnetwork(n, k, q) 
		# net, E = scalefreenetwork(n, 2*n, 4)
		sources, sinks = sourcesinklocs(ns, nd, n)

		#Make the initial flow 
		m = numedges(net)
		X = zeros(Float64, m) .+ 1.0
		P = sourcesinkvector(sources, sinks, n, Ptot)
		L = weightedlaplacianmatrix(E, X) 
		theta = nodevoltages(L, P)
		F = edgecurrents(theta, E, X)
		flowmaxtemp = maximum(abs.(F))
		F = F./flowmaxtemp 

		for j=1:subensemblesize
			ind = (i-1)*subensemblesize + j
			ensemblemember = CascadeNetwork(net, E, sources, sinks, F) 
			ensemble[ind] = ensemblemember 
			triggervec[ind] = triggerindex(F)
			degs = getdegrees(E)
			degmeans[ind] = mean(degs) 
			fmaxes[ind] = flowmaxtemp 
		end
	end

	alpha = 0.01
	stepsize = 0.3
	beneath0 = true 
	beneath1 = true 
	iterations = 1 
	ensemblesize = supensemblesize*subensemblesize 
	while abs(stepsize) > tol 
		println(alpha)

		#For every member in the ensemble
		runningsurvivors = Array{Float64, 1}(undef, ensemblesize) 
		for i=1:ensemblesize
			# println(i)


			#Unpack the network 
			net = ensemble[i].network 
			E = ensemble[i].E 
			m = numedges(net)
			X = zeros(Float64, m) .+ 1.0
			sources = ensemble[i].sources 
			sinks = ensemble[i].sinks
			F = ensemble[i].initialflow 
			nodeset = collect(1:n)
			edgeset = collect(1:m)
			P = sourcesinkvector(sources, sinks, n, Ptot)

			#Reactivate all edges 
			for j=1:length(net.edgelist)
				net.edgelist[j].active = true 
			end

			#Knock out trigger
			d = triggervec[i]
			edgeset1 = filter(x->x≠d, edgeset) 
			E1 = E[1:end, 1:end .≠ d]
			net.edgelist[d].active = false 

			#Check to see if its broken 
			Adj = adjacencymatrix(E1)
			components, table = connectedcomponents(Adj)

			#Do cascade 
			flowvec = zeros(Float64, length(F)) .+ F   
			flowvec[d] = NaN
			info = newrecorder(m, ns, Ptot, nd, sum(abs.(P))/2.0, flowvec) 
			simtime = 0
			for j=1:components
				chunk = table[j] 
				nodesubset = zeros(Int64, length(chunk))
				for l=1:length(chunk)
					nodesubset[l] = nodeset[chunk[l]]
				end
				info = fracture!(net, nodesubset, P, X, F, info, simtime, alpha, fmaxes[i], false)
			end  

			#Record fraction of edges remaining
			T, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, flowseq = cascaderesults(info) 
			runningsurvivors[i] = fracedges

		end 

		#If fraction has crossed the target then reverse and half the step size 
		global avremaining = mean(runningsurvivors) 
		if avremaining > 0.8
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

	return alpha, mean(fmaxes), iterations, degmeans  

end 



function sequencedstartbisection(tol::Float64, n::Int64, ns::Int64, nd::Int64, q::Float64, k::Int64) 


	#Construct a network
	net, E = smallworldnetwork(n, k, q) 
	# net, E = scalefreenetwork(n, 2*n, 3)
	sources, sinks = sourcesinklocs(ns, nd, n)
	# println("done")

	#Make the initial flow 
	m = numedges(net)
	# println(m)
	X = zeros(Float64, m) .+ 1.0
	P = sourcesinkvector(sources, sinks, n, 1.0)
	L = weightedlaplacianmatrix(E, X)  
	theta = nodevoltages(L, P)
	F = edgecurrents(theta, E, X)
	fmax = maximum(abs.(F))
	F = F./fmax 
	nodeset = collect(1:n)
	edgeset = collect(1:m)  

	alphacs = zeros(Float64, m) 
	relativeload = zeros(Float64, m) 

	#Iterate through triggering edges 
	for z=1:m 

		#Start bisection 
		alpha = 0.01
		stepsize = 0.3 
		beneath0 = true 
		beneath1 = true 
		iterations = 1 
		
		relativeload[z] = abs(F[z])/mean(abs.(F)) 
		
		while abs(stepsize) > tol 

			#Reactivate all edges 
			for i=1:m
				net.edgelist[i].active = true 
			end

			#Knock out trigger edge  
			edgeset1 = filter(x->x≠z, edgeset) 
			E1 = E[1:end, 1:end .≠ z]
			net.edgelist[z].active = false 

			#Check to see if its broken 
			Adj = adjacencymatrix(E1)
			components, table = connectedcomponents(Adj)

			#Do cascade 
			flowvec = zeros(Float64, length(F)) .+ F   
			flowvec[z] = NaN
			info = newrecorder(m, ns, 1.0, nd, sum(abs.(P))/2.0, flowvec) 
			simtime = 0
			for j=1:components
				chunk = table[j] 
				nodesubset = zeros(Int64, length(chunk))
				for l=1:length(chunk)
					nodesubset[l] = nodeset[chunk[l]]
				end
				info = fracture!(net, nodesubset, P, X, F, info, simtime, alpha, fmax, false)
			end  

			#Record fraction of edges remaining
			T, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, flowseq = cascaderesults(info) 
			
			#If fraction has crossed the target then reverse and half the step size 
			if fracedges > 0.5
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
		alphacs[z] = alpha 
	end 
	return relativeload, alphacs, fmax 
end 



function cascadebisection(ensemblesize::Int64, tol::Float64, n::Int64, ns::Int64, nd::Int64, q::Float64, k::Int64, Ptot::Float64) 


	#Make an ensemble of networks here, with node placements
	#and with initial flow pattern
	ensemble = Array{CascadeNetwork, 1}(undef, ensemblesize)
	fmaxes = Array{Float64, 1}(undef, ensemblesize) 
	degmeans = Array{Float64, 1}(undef, ensemblesize) 
	for i=1:ensemblesize

		#Construct each network
		net, E = smallworldnetwork(n, k, q) 
		# net, E = scalefreenetwork(n, 2*n, 4)
		degs = getdegrees(E)
		degmeans[i] = mean(degs) 
		sources, sinks = sourcesinklocs(ns, nd, n)

		#Make the initial flow 
		m = numedges(net)
		X = zeros(Float64, m) .+ 1.0
		P = sourcesinkvector(sources, sinks, n, Ptot)
		L = weightedlaplacianmatrix(E, X) 
		theta = nodevoltages(L, P)
		F = edgecurrents(theta, E, X)
		flowmaxtemp = maximum(abs.(F))
		fmaxes[i] = flowmaxtemp 
		F = F./flowmaxtemp 

		#Store the ensemble member 
		ensemblemember = CascadeNetwork(net, E, sources, sinks, F) 
		ensemble[i] = ensemblemember 
		# println(i) 
	end

	alpha = 0.01
	stepsize = 0.3
	beneath0 = true 
	beneath1 = true 
	iterations = 1 
	while abs(stepsize) > tol 
		# println(alpha)

		#For every member in the ensemble
		runningsurvivors = Array{Float64, 1}(undef, ensemblesize) 
		for i=1:ensemblesize


			#Unpack the network 
			net = ensemble[i].network 
			E = ensemble[i].E 
			m = numedges(net)
			X = zeros(Float64, m) .+ 1.0
			sources = ensemble[i].sources 
			sinks = ensemble[i].sinks
			F = ensemble[i].initialflow 
			nodeset = collect(1:n)
			edgeset = collect(1:m)
			P = sourcesinkvector(sources, sinks, n, Ptot)

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
			flowvec = zeros(Float64, length(F)) .+ F   
			flowvec[d] = NaN
			info = newrecorder(m, ns, Ptot, nd, sum(abs.(P))/2.0, flowvec) 
			simtime = 0
			for j=1:components
				chunk = table[j] 
				nodesubset = zeros(Int64, length(chunk))
				for l=1:length(chunk)
					nodesubset[l] = nodeset[chunk[l]]
				end
				info = fracture!(net, nodesubset, P, X, F, info, simtime, alpha, fmaxes[i], false)
			end  

			#Record fraction of edges remaining
			T, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, flowseq = cascaderesults(info) 
			runningsurvivors[i] = fracedges

		end 

		#If fraction has crossed the target then reverse and half the step size 
		global avremaining = mean(runningsurvivors) 
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

	return alpha, mean(fmaxes), iterations, degmeans  

end 


function cascadebisection2(ensemblesize::Int64, tol::Float64, n::Int64, ns::Int64, nd::Int64, q::Float64, k::Int64, Ptot::Float64) 


	#Make an ensemble of networks here, with node placements
	#and with initial flow pattern
	ensemble = Array{CascadeNetwork, 1}(undef, ensemblesize)
	fmaxes = Array{Float64, 1}(undef, ensemblesize) 
	degmeans = Array{Float64, 1}(undef, ensemblesize) 
	for i=1:ensemblesize

		#Construct each network
		net, E = smallworldnetwork(n, k, q) 
		# net, E = scalefreenetwork(n, 2*n, 3)
		degs = getdegrees(E)
		degmeans[i] = mean(degs) 
		sources, sinks = sourcesinklocs(ns, nd, n)

		#Make the initial flow 
		m = numedges(net)
		X = zeros(Float64, m) .+ 1.0
		P = sourcesinkvector(sources, sinks, n, Ptot)
		L = weightedlaplacianmatrix(E, X) 
		theta = nodevoltages(L, P)
		F = edgecurrents(theta, E, X)
		flowmaxtemp = maximum(abs.(F))
		fmaxes[i] = flowmaxtemp 
		F = F./flowmaxtemp 

		#Store the ensemble member 
		ensemblemember = CascadeNetwork(net, E, sources, sinks, F) 
		ensemble[i] = ensemblemember 
		# println(i) 
	end

	alpha = 0.01
	stepsize = 0.3
	beneath0 = true 
	beneath1 = true 
	iterations = 1 
	while abs(stepsize) > tol 
		# println(alpha)

		#For every member in the ensemble
		runningsurvivors = Array{Float64, 1}(undef, ensemblesize) 
		for i=1:ensemblesize


			#Unpack the network 
			net = ensemble[i].network 
			E = ensemble[i].E 
			m = numedges(net)
			X = zeros(Float64, m) .+ 1.0
			sources = ensemble[i].sources 
			sinks = ensemble[i].sinks
			F = ensemble[i].initialflow 
			nodeset = collect(1:n)
			edgeset = collect(1:m)
			P = sourcesinkvector(sources, sinks, n, Ptot)

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
			flowvec = zeros(Float64, length(F)) .+ F   
			flowvec[d] = NaN
			info = newrecorder(m, ns, Ptot, nd, sum(abs.(P))/2.0, flowvec) 
			simtime = 0
			remainingedges = 0 
			for j=1:components
				chunk = table[j] 
				nodesubset = zeros(Int64, length(chunk))
				for l=1:length(chunk)
					nodesubset[l] = nodeset[chunk[l]]
				end
				remainingedges += fracture2!(net, nodesubset, P, X, F, simtime, alpha, fmaxes[i])
			end  

			#Record fraction of edges remaining
			runningsurvivors[i] = remainingedges/Float64(m) 

		end 

		#If fraction has crossed the target then reverse and half the step size 
		global avremaining = mean(runningsurvivors) 
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

	return alpha, mean(fmaxes), iterations, degmeans  
end 


function cascadebisection3rand(ensemblesize::Int64, tol::Float64, n::Int64, ns::Int64, nd::Int64, q::Float64, k::Int64, Ptot::Float64) 


	#Make an ensemble of networks here, with node placements
	#and with initial flow pattern
	ensemble = Array{CascadeNetworkRand, 1}(undef, ensemblesize)
	fmaxes = Array{Float64, 1}(undef, ensemblesize)
	degmeans = Array{Float64, 1}(undef, ensemblesize) 
	for i=1:ensemblesize

		#Construct each network
		net, E = smallworldnetwork(n, k, q) 
		# net, E = scalefreenetwork(n, 2*n, 3)
		degs = getdegrees(E)
		degmeans[i] = mean(degs) 
		sources, sinks = sourcesinklocs(ns, nd, n)

		#Make the initial flow 
		m = numedges(net)
		X = zeros(Float64, m) .+ 1.0
		P = sourcesinkvector(sources, sinks, n, Ptot)
		L = weightedlaplacianmatrix(E, X) 
		theta = nodevoltages(L, P)
		F = edgecurrents(theta, E, X)
		flowmaxtemp = maximum(abs.(F))
		fmaxes[i] = flowmaxtemp 
		F = F./flowmaxtemp 

		#Make noise for the capacities 
		d = Normal(0.0, 0.05)
		alphanoise = rand(d, m) 

		#Store the ensemble member 
		ensemblemember = CascadeNetworkRand(net, E, sources, sinks, F, alphanoise) 
		ensemble[i] = ensemblemember 
		# println(i) 
	end

	alpha = 0.01
	stepsize = 0.3
	beneath0 = true 
	beneath1 = true 
	iterations = 1 
	while abs(stepsize) > tol 
		# println(alpha)

		#For every member in the ensemble
		runningsurvivors = Array{Float64, 1}(undef, ensemblesize) 
		for i=1:ensemblesize


			#Unpack the network 
			net = ensemble[i].network 
			E = ensemble[i].E 
			capnoise = ensemble[i].alphanoise
			m = numedges(net)
			X = zeros(Float64, m) .+ 1.0
			sources = ensemble[i].sources 
			sinks = ensemble[i].sinks
			F = ensemble[i].initialflow 
			nodeset = collect(1:n)
			edgeset = collect(1:m)
			P = sourcesinkvector(sources, sinks, n, Ptot)

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
			flowvec = zeros(Float64, length(F)) .+ F   
			flowvec[d] = NaN
			info = newrecorder(m, ns, Ptot, nd, sum(abs.(P))/2.0, flowvec) 
			simtime = 0
			remainingedges = 0 
			for j=1:components
				chunk = table[j] 
				nodesubset = zeros(Int64, length(chunk))
				for l=1:length(chunk)
					nodesubset[l] = nodeset[chunk[l]]
				end
				remainingedges += fracturerand!(net, nodesubset, P, X, F, capnoise, simtime, alpha, fmaxes[i])
			end  

			#Record fraction of edges remaining
			runningsurvivors[i] = remainingedges/Float64(m) 

		end 

		#If fraction has crossed the target then reverse and half the step size 
		global avremaining = mean(runningsurvivors) 
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

	return alpha, mean(fmaxes), iterations, degmeans  
end 


function smcascadefixedtop(ensemblesize::Int64, n::Int64, ns::Int64, nd::Int64, q::Float64, alphares::Int64, alphamin::Float64, alphamax::Float64, k::Int64, Ptot::Float64, Xval::Float64)
	#A single network with an ensemble of different source-sink locations

	#Create a network topology
	net, E0 =  austriannetwork() #smallworldnetwork(n, k, q) 
	m = numedges(net)
	nodeset = collect(1:n)
	edgeset = collect(1:m)  
	X = zeros(Float64, m) .+ Xval 

	alphas = collect(range(alphamin, length=alphares, stop=alphamax))
	edgemat = zeros(Float64, ensemblesize, alphares)  
	Tmat = zeros(Float64, ensemblesize, alphares) 

	#Loop over ensemble
	flowmax = zeros(Float64, ensemblesize)
	for z=1:ensemblesize
		# println(z)
		
		#Create a P vector
		sources, sinks = sourcesinklocsmod(net, ns, nd, n)
		P = noisysourcesinkvector(sources, sinks, n, Ptot)

		#Solve the DC flow equations
		L = weightedlaplacianmatrix(E0, 1.0./X)  
		theta = nodevoltages(L, P)
		F = edgecurrents(theta, E0, X)
		flowmaxtemp = maximum(abs.(F))
		F = F./flowmaxtemp 
		flowmax[z] = flowmaxtemp

		#Empty vectors to record data 
		rT = zeros(Float64, alphares)
		redges = zeros(Float64, alphares)
		rsources = zeros(Float64, alphares)
		rsinks = zeros(Float64, alphares)
		rweightedsources = zeros(Float64, alphares)
		rpower = zeros(Float64, alphares)

		#Loop over alpha 
		# println(z)
		for (ind, alpha) ∈ enumerate(alphas)
			# println(ind)

			#Reactive all the edges
			for i=1:length(net.edgelist)
				net.edgelist[i].active = true 
			end

			#Knock out most loaded edge 
			d = argmax(abs.(F)) 
			edgeset1 = filter(x->x≠d, edgeset) 
			E1 = E0[1:end, 1:end .≠ d]
			net.edgelist[d].active = false 

			#Check to see if its broken 
			Adj = adjacencymatrix(E1)
			components, table = connectedcomponents(Adj)

			flowvec = zeros(Float64, length(F)) .+ F   
			flowvec[d] = NaN
			info = newrecorder(m, ns, Ptot, nd, sum(abs.(P))/2.0, flowvec) 
			simtime = 0 

			for i=1:components
				chunk = table[i] 
				nodesubset = zeros(Int64, length(chunk))
				for j=1:length(chunk)
					nodesubset[j] = nodeset[chunk[j]]
				end
				info = fracture!(net, nodesubset, P, X, F, info, simtime, alpha, flowmaxtemp, false)
			end

			#Build up vectors of rT etc
			T, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, flowseq = cascaderesults(info)
			rT[ind] = T
			redges[ind] = fracedges
			rsources[ind] = fracsources
			rsinks[ind] = fracsinks
			rweightedsources[ind] = fracweightedsources
			rpower[ind] = fracpower
		end

		#Stack up vecs of rT etc.
		edgemat[z,:] = redges
		Tmat[z,:] = rpower
	end 

	# p0 = 10
	# p1 = 15 
	# p2 = 75

	#return the mean of the rT etc data 
	return alphas, mean(edgemat, dims=1)[1,:], std(edgemat, dims=1)[1,:], mean(Tmat, dims=1)[1,:], std(Tmat, dims=1)[1,:], edgemat, Tmat, mean(flowmax)
end


function simplecascade(n::Int64, ns::Int64, nd::Int64, q::Float64, alpha::Float64, k::Int64, Ptot::Float64, Xval::Float64)

	#Create a network topology
	net, E0 = smallworldnetwork(n, k, q) 
	m = numedges(net)
	nodeset = collect(1:n)
	edgeset = collect(1:m)  
	X = zeros(Float64, m) .+ Xval

	#Create a P vector
	sources, sinks = sourcesinklocs(ns, nd, n)
	P = sourcesinkvector(sources, sinks, n, Ptot*n)

	#Solve the DC flow equations
	L = weightedlaplacianmatrix(E0, 1.0./X)  
	theta = nodevoltages(L, P)
	F = edgecurrents(theta, E0, X) 
	flowmaxtemp = maximum(abs.(F))
	F = F./flowmaxtemp 
	meanflow = mean(abs.(F))

	#Reactive all the edges
	for i=1:length(net.edgelist)
		net.edgelist[i].active = true 
	end

	d = argmax(abs.(F)) 
	edgeset1 = filter(x->x≠d, edgeset) 
	E1 = E0[1:end, 1:end .≠ d]
	net.edgelist[d].active = false 

	#Check to see if its broken 
	Adj = adjacencymatrix(E1)
	components, table = connectedcomponents(Adj)

	fragments = 0
	for i=1:components
		chunk = table[i] 
		nodesubset = zeros(Int64, length(chunk))
		for j=1:length(chunk)
			nodesubset[j] = nodeset[chunk[j]]
		end
		fragments = fracture3!(net, nodesubset, P, X, F, 0, alpha, flowmaxtemp)
		totalfragments += fragments
	end

	return totalfragments
end


function smcascadedist(ensemblesize::Int64, n::Int64, ns::Int64, nd::Int64, q::Float64, alpha::Float64, k::Int64, Ptot::Float64, Xval::Float64)
	#A single network with an ensemble of different source-sink locations

	#Create a network topology
	net, E0 = smallworldnetwork(n, k, q) 
	m = numedges(net)
	nodeset = collect(1:n)
	edgeset = collect(1:m)  
	X = zeros(Float64, m) .+ Xval 

	#Create a P vector
	sources, sinks = sourcesinklocs(ns, nd, n)
	P = sourcesinkvector(sources, sinks, n, Ptot*n)

	#Solve the DC flow equations
	L = weightedlaplacianmatrix(E0, 1.0./X)  
	theta = nodevoltages(L, P)
	F = edgecurrents(theta, E0, X) 
	flowmaxtemp = maximum(abs.(F))
	F = F./flowmaxtemp 
	meanflow = mean(abs.(F))

	relativeload = zeros(Float64, ensemblesize)
	survivingedges = zeros(Float64, ensemblesize) 
	survivingpower = zeros(Float64, ensemblesize) 

	#Loop over ensemble
	for z=1:ensemblesize
		# println(z)

		#Reactive all the edges
		for i=1:length(net.edgelist)
			net.edgelist[i].active = true 
		end

		#Knock out most loaded edge 
		d = rand(1:m)
		# d = argmax(abs.(F)) 
		edgeset1 = filter(x->x≠d, edgeset) 
		E1 = E0[1:end, 1:end .≠ d]
		net.edgelist[d].active = false 

		#Compute and save relative load of knocked out edge 
		relativeload[z] = abs(F[d])/meanflow
		
		#Check to see if its broken 
		Adj = adjacencymatrix(E1)
		components, table = connectedcomponents(Adj)

		flowvec = zeros(Float64, length(F)) .+ F   
		flowvec[d] = NaN
		info = newrecorder(m, ns, Ptot, nd, sum(abs.(P))/2.0, flowvec) 
		simtime = 0 

		for i=1:components
			chunk = table[i] 
			nodesubset = zeros(Int64, length(chunk))
			for j=1:length(chunk)
				nodesubset[j] = nodeset[chunk[j]]
			end
			info = fracture!(net, nodesubset, P, X, F, info, simtime, alpha, flowmaxtemp, false)
		end

		#Build up vectors of rT etc
		T, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, flowseq = cascaderesults(info)

		survivingedges[z] = fracedges 
		survivingpower[z] = fracpower


	end 

	#return the mean of the rT etc data 
	return relativeload, survivingedges, survivingpower 
end



function smcascadeclock(n::Int64, pos::Int64, q::Float64, alphares::Int64, alphamin::Float64, alphamax::Float64, k::Int64, Ptot::Float64, Xval::Float64)
	#A single network with with a fixed source-sink pair location
	#Create a network topology

	ensemblesize = 1
	net, E0 = smallworldnetwork(n, k, q) 
	m = numedges(net)
	nodeset = collect(1:n)
	edgeset = collect(1:m)  
	X = zeros(Float64, m) .+ Xval 

	alphas = collect(range(alphamin, length=alphares, stop=alphamax))
	edgemat = zeros(Float64, ensemblesize, alphares)  
	Tmat = zeros(Float64, ensemblesize, alphares) 

	#Loop over ensemble
	for z=1:ensemblesize
		
		#Create a P vector
		sources, sinks = sourcessinkclockface(n, pos)
		P = sourcesinkvector(sources, sinks, n, Ptot)

		#Solve the DC flow equations
		L = weightedlaplacianmatrix(E0, 1.0./X)  
		theta = nodevoltages(L, P)
		F = edgecurrents(theta, E0, X) 
		flowmaxtemp = maximum(abs.(F))
		F = F./flowmaxtemp 
		meanflow = mean(abs.(F))

		#Empty vectors to record data 
		rT = zeros(Float64, alphares)
		redges = zeros(Float64, alphares)
		rsources = zeros(Float64, alphares)
		rsinks = zeros(Float64, alphares)
		rweightedsources = zeros(Float64, alphares)
		rpower = zeros(Float64, alphares)

		#Loop over alpha 
		for (ind, alpha) ∈ enumerate(alphas)

			#Reactive all the edges
			for i=1:length(net.edgelist)
				net.edgelist[i].active = true 
			end

			#Knock out most loaded edge 
			d = argmax(abs.(F)) 
			edgeset1 = filter(x->x≠d, edgeset) 
			E1 = E0[1:end, 1:end .≠ d]
			net.edgelist[d].active = false 

			#Check to see if its broken 
			Adj = adjacencymatrix(E1)
			components, table = connectedcomponents(Adj)

			flowvec = zeros(Float64, length(F)) .+ abs.(F) 
			flowvec[d] = NaN
			info = newrecorder(m, 2, Ptot, 62, sum(abs.(P))/2.0, flowvec) 
			simtime = 0 

			for i=1:components
				chunk = table[i] 
				nodesubset = zeros(Int64, length(chunk))
				for j=1:length(chunk)
					nodesubset[j] = nodeset[chunk[j]]
				end
				info = fracture!(net, nodesubset, P, X, F, info, simtime, alpha, flowmaxtemp, false)
			end

			#Build up vectors of rT etc
			T, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, flowseq = cascaderesults(info)
			rT[ind] = T
			redges[ind] = fracedges
			rsources[ind] = fracsources
			rsinks[ind] = fracsinks
			rweightedsources[ind] = fracweightedsources
			rpower[ind] = fracpower
		end

		#Stack up vecs of rT etc.
		edgemat[z,:] = redges
		Tmat[z,:] = rT
	end 

	#return the mean of the rT etc data 
	return alphas, mean(edgemat, dims=1)[1,:], std(edgemat, dims=1)[1,:], mean(Tmat, dims=1)[1,:], std(Tmat, dims=1)[1,:], edgemat
end


function smcascadeclock2(n::Int64, pos::Int64, q::Float64, alphares::Int64, alphamin::Float64, alphamax::Float64, k::Int64, Ptot::Float64, Xval::Float64)
	#A single network with with a fixed source-sink pair location
	#Create a network topology

	ensemblesize = 1
	net, E0 = smallworldnetwork(n, k, q) 
	m = numedges(net)
	nodeset = collect(1:n)
	edgeset = collect(1:m)  
	X = zeros(Float64, m) .+ Xval 

	alphas = collect(range(alphamin, length=alphares, stop=alphamax))
	edgemat = zeros(Float64, alphares)  

		
	#Create a P vector
	sources, sinks = sourcessinkclockface(n, pos)
	P = sourcesinkvector(sources, sinks, n, Ptot)

	#Solve the DC flow equations
	L = weightedlaplacianmatrix(E0, 1.0./X)  
	theta = nodevoltages(L, P)
	F = edgecurrents(theta, E0, X) 
	flowmaxtemp = maximum(abs.(F))
	F = F./flowmaxtemp 
	meanflow = mean(abs.(F))

	#Loop over alpha 
	for (ind, alpha) ∈ enumerate(alphas)

		#Reactive all the edges
		for i=1:length(net.edgelist)
			net.edgelist[i].active = true 
		end

		#Knock out most loaded edge 
		d = argmax(abs.(F)) 
		edgeset1 = filter(x->x≠d, edgeset) 
		E1 = E0[1:end, 1:end .≠ d]
		net.edgelist[d].active = false 

		#Check to see if its broken 
		Adj = adjacencymatrix(E1)
		components, table = connectedcomponents(Adj)

		flowvec = zeros(Float64, length(F)) .+ abs.(F) 
		simtime = 0 

		for i=1:components
			chunk = table[i] 
			nodesubset = zeros(Int64, length(chunk))
			for j=1:length(chunk)
				nodesubset[j] = nodeset[chunk[j]]
			end
			edgemat[ind] = fracture2!(net, nodesubset, P, X, F, simtime, alpha, flowmaxtemp)
		end
	end 

	#return the mean of the rT etc data 
	return edgemat 
end


function smcascadeclockrand(n::Int64, pos::Int64, q::Float64, alphares::Int64, alphamin::Float64, alphamax::Float64, k::Int64, Ptot::Float64, Xval::Float64)
	#A single network with with a fixed source-sink pair location
	#Create a network topology

	ensemblesize = 1
	net, E0 = smallworldnetwork(n, k, q) 
	m = numedges(net)
	nodeset = collect(1:n)
	edgeset = collect(1:m)  
	X = zeros(Float64, m) .+ Xval 

	alphas = collect(range(alphamin, length=alphares, stop=alphamax))
	edgemat = zeros(Float64, alphares)  

	#Make edge noise 
	d = Normal(0.0, 0.05)
	alphanoise = rand(d, m)
		
	#Create a P vector
	sources, sinks = sourcessinkclockface(n, pos)
	P = sourcesinkvector(sources, sinks, n, Ptot)

	#Solve the DC flow equations
	L = weightedlaplacianmatrix(E0, 1.0./X)  
	theta = nodevoltages(L, P)
	F = edgecurrents(theta, E0, X) 
	flowmaxtemp = maximum(abs.(F))
	F = F./flowmaxtemp 
	meanflow = mean(abs.(F))

	#Loop over alpha 
	for (ind, alpha) ∈ enumerate(alphas)

		#Reactive all the edges
		for i=1:length(net.edgelist)
			net.edgelist[i].active = true 
		end

		#Knock out most loaded edge 
		d = argmax(abs.(F)) 
		edgeset1 = filter(x->x≠d, edgeset) 
		E1 = E0[1:end, 1:end .≠ d]
		net.edgelist[d].active = false 

		#Check to see if its broken 
		Adj = adjacencymatrix(E1)
		components, table = connectedcomponents(Adj)

		flowvec = zeros(Float64, length(F)) .+ abs.(F) 
		simtime = 0 

		for i=1:components
			chunk = table[i] 
			nodesubset = zeros(Int64, length(chunk))
			for j=1:length(chunk)
				nodesubset[j] = nodeset[chunk[j]]
			end
			edgemat[ind] = fracturerand!(net, nodesubset, P, X, F, alphanoise, simtime, alpha, flowmaxtemp)
		end
	end 

	#return the mean of the rT etc data 
	return edgemat 
end


# function cascadesequence(net::Graph, E::Matrix{Float64}, sources::Vector{Int64}, sinks::Vector{Int64}, alpha::Float64, Ptot::Float64, Xval::Float64)

# 	n = size(E,1)
# 	m = size(E,2) 
# 	nodeset = collect(1:n)
# 	edgeset = collect(1:m)
# 	X = zeros(Float64, m) .+ Xval
# 	ns = length(sources)
# 	nd = length(sinks) 
# 	P = sourcesinkvector(sources, sinks, n, Ptot*n)

# 	#Solve the DC flow equations
# 	L = weightedlaplacianmatrix(E, 1.0./X)  
# 	theta = nodevoltages(L, P)
# 	F = edgecurrents(theta, E, X)

# 	#Active all the edges
# 	for i=1:m
# 		net.edgelist[i].active = true 
# 	end

# 	#Knock out most loaded edge 
# 	d = argmax(abs.(F)) 
# 	edgeset1 = filter(x->x≠d, edgeset) 
# 	E1 = E[1:end, 1:end .≠ d]
# 	net.edgelist[d].active = false 

# 	#Check to see if its broken 
# 	Adj = adjacencymatrix(E1)
# 	components, table = connectedcomponents(Adj)

# 	flowvec = zeros(Float64, length(F)) .+ abs.(F)  
# 	println("Initial flow:")
# 	println(flowvec)
# 	flowvec[d] = NaN
# 	info = newrecorder(m, ns, Ptot, nd, sum(abs.(P))/2.0, flowvec) 
# 	simtime = 0 

# 	for i=1:components
# 		chunk = table[i] 
# 		nodesubset = zeros(Int64, length(chunk))
# 		for j=1:length(chunk)
# 			nodesubset[j] = nodeset[chunk[j]]
# 		end
# 		info = fracture!(net, nodesubset, P, X, F, info, simtime, alpha, true)
# 	end

# 	T, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, flowseq = cascaderesults(info)
# 	return T, fracedges, flowseq
# end


# function singlecascade(net::Graph, E0::Matrix{Float64}, ensemblesize::Int64, n::Int64, ns::Int64, nd::Int64, q::Float64, alphares::Int64, alphamin::Float64, alphamax::Float64, k::Int64, Ptot::Float64, Xval::Float64)
# 	#A single network with an ensemble of different source-sink locations

# 	#Create a network topology
# 	m = numedges(net)
# 	nodeset = collect(1:n)
# 	edgeset = collect(1:m)  
# 	X = zeros(Float64, m) .+ Xval 

# 	alphas = collect(range(alphamin, length=alphares, stop=alphamax))
# 	edgemat = zeros(Float64, ensemblesize, alphares)  
# 	Tmat = zeros(Float64, ensemblesize, alphares) 

# 	#Loop over ensemble
# 	flowmax = zeros(Float64, ensemblesize)
# 	for z=1:ensemblesize
# 		# println(z)
		
# 		#Create a P vector
# 		sources, sinks = sourcesinklocs(ns, nd, n)
# 		P = sourcesinkvector(sources, sinks, n, Ptot*n)

# 		#Solve the DC flow equations
# 		L = weightedlaplacianmatrix(E0, 1.0./X)  
# 		theta = nodevoltages(L, P)
# 		F = edgecurrents(theta, E0, X)
# 		flowmax[z] = maximum(abs.(F))

# 		#Empty vectors to record data 
# 		rT = zeros(Float64, alphares)
# 		redges = zeros(Float64, alphares)
# 		rsources = zeros(Float64, alphares)
# 		rsinks = zeros(Float64, alphares)
# 		rweightedsources = zeros(Float64, alphares)
# 		rpower = zeros(Float64, alphares)

# 		#Loop over alpha 
# 		# println(z)
# 		for (ind, alpha) ∈ enumerate(alphas)
# 			# println(ind)

# 			#Reactive all the edges
# 			for i=1:length(net.edgelist)
# 				net.edgelist[i].active = true 
# 			end

# 			#Knock out most loaded edge 
# 			d = argmax(abs.(F)) 
# 			edgeset1 = filter(x->x≠d, edgeset) 
# 			E1 = E0[1:end, 1:end .≠ d]
# 			net.edgelist[d].active = false 

# 			#Check to see if its broken 
# 			Adj = adjacencymatrix(E1)
# 			components, table = connectedcomponents(Adj)

# 			flowvec = zeros(Float64, length(F)) .+ F   
# 			flowvec[d] = NaN
# 			info = newrecorder(m, ns, Ptot, nd, sum(abs.(P))/2.0, flowvec) 
# 			simtime = 0 

# 			for i=1:components
# 				chunk = table[i] 
# 				nodesubset = zeros(Int64, length(chunk))
# 				for j=1:length(chunk)
# 					nodesubset[j] = nodeset[chunk[j]]
# 				end
# 				info = fracture!(net, nodesubset, P, X, F, info, simtime, alpha, false)
# 			end

# 			#Build up vectors of rT etc
# 			T, fracedges, fracsources, fracweightedsources, fracsinks, fracpower, flowseq = cascaderesults(info)
# 			rT[ind] = T
# 			redges[ind] = fracedges
# 			rsources[ind] = fracsources
# 			rsinks[ind] = fracsinks
# 			rweightedsources[ind] = fracweightedsources
# 			rpower[ind] = fracpower
# 		end

# 		#Stack up vecs of rT etc.
# 		edgemat[z,:] = redges
# 		Tmat[z,:] = rT
# 	end 

# 	# p0 = 10
# 	# p1 = 15 
# 	# p2 = 75

# 	#return the mean of the rT etc data 
# 	return alphas, mean(edgemat, dims=1)[1,:], std(edgemat, dims=1)[1,:], mean(Tmat, dims=1)[1,:], std(Tmat, dims=1)[1,:], edgemat, Tmat, mean(flowmax)
# end














































