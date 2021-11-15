include("graphs.jl") 
using DelimitedFiles

function incidencematrix(graph::Graph)
	n = numnodes(graph)
	m = numedges(graph)
	E = zeros(Float64, n, m)
	for i=1:n
		v = graph.nodelist[i]
		for j in v.outedges
			E[i,j] = 1.0 
		end
		for j in v.inedges
			E[i,j] = -1.0 
		end
	end
	return E 
end


function incidencematrixtonetwork(E::Matrix{Float64})
	n = size(E,1)
	m = size(E,2) 
	net = Graph(Node[], Edge[], 1, 1)
	for i=1:n 
		addnode!(net)
	end
	for j=1:m 
		origin = 0 
		target = 0 
		for i=1:n 
			if E[i,j] < -0.1 
				origin = i
			end 
			if E[i,j] > 0.1
				target = i 
			end
		end
		addedge!(net, origin, target) 
	end
	return net 
end


function getdegrees(E::Matrix{Float64})
	n = size(E,1) 
	m = size(E,2)  
	degrees = zeros(Float64, n)
	for i=1:n 
		for j=1:m
			degrees[i] += abs(E[i,j]) 
		end
	end
	return degrees 
end


function adjacencymatrix(E::Matrix{Float64})
	return Matrix(Diagonal(getdegrees(E))) - laplacianmatrix(E) 
end


function laplacianmatrix(E::Matrix{Float64})
	return E*transpose(E)
end


function weightedlaplacianmatrix(E::Matrix{Float64}, w::Vector{Float64})
	return E*Matrix(Diagonal(w))*transpose(E) 
end


function wattsstrogatz!(graph::Graph, E::Matrix{Float64}, q::Float64) 
	n = size(E,1) 
	m = size(E,2)
	A = adjacencymatrix(E) 
	for j=1:m
		z = rand() 
		if z < q
			e = graph.edgelist[j]
			s = e.origin
			t = e.target
			while true 
				possible_t = rand(1:n)
				if (possible_t ≠ s) && (A[s, possible_t] < 0.1)
					A[s, possible_t] = 1.0
					A[possible_t, s] = 1.0
					A[s, t] = 0.0
					A[t, s] = 0.0
					E[t, j] = 0.0
					E[possible_t, j] = -1.0
					break 
				end 
			end
		end
	end
end


function scalefreenetwork(n::Int64, m::Int64, m0::Int64)
	net = Graph(Node[], Edge[], 1, 1)
	for i=1:4
		addnode!(net)
	end
	addedge!(net, 1, 2)
	addedge!(net, 1, 3)
	addedge!(net, 3, 4)
	addedge!(net, 1, 4) 

	while numnodes(net) < n 
		addnode!(net) 
		v = net.nodelist[end] 
		i = v.nodeindex 
		while length(v.outedges) < m0 
			j = rand(1:i-1)
			w = net.nodelist[j]
			k = length(w.outedges) + length(w.inedges) 
			if rand() < Float64(k)/(2.0*numedges(net))
				addedge!(net, i, j)
			end 
		end
	end

	E = incidencematrix(net)
	return net, E 
end


function smallworldnetwork(n::Int64, k::Int64, q::Float64)
	net = Graph(Node[], Edge[], 1, 1)
	for i=1:n 
		addnode!(net) 
	end
	if k ≥ n-1
		println("complete")
		for i=1:(n-1)  
			for j=(i+1):n 
				addedge!(net, i, j) 
			end
		end
	else
		for i=1:n
			for j=1:div(k,2)
				addedge!(net, i, mod1((i+j),n)) 
			end
		end
	end 
	E = incidencematrix(net)
	wattsstrogatz!(net, E, q) 
	net = incidencematrixtonetwork(E)  
	return net, E 
end 


function austriannetwork()
	E = readdlm("austria_incidence.txt")
	net = incidencematrixtonetwork(E) 
	return net, E 
end


function testnetwork()
	net = Graph(Node[], Edge[], 1, 1)
	for i=1:6
		addnode!(net) 
	end
	addedge!(net, 1, 2)
	addedge!(net, 1, 3)
	addedge!(net, 2, 3)
	addedge!(net, 3, 4)
	addedge!(net, 4, 5)
	addedge!(net, 4, 6)
	addedge!(net, 5, 6)
	E = incidencematrix(net) 
	X = zeros(Float64, 7)
	X .+= 1.0 
	P = zeros(Float64, 6) 
	P[1] = 1.0 
	P[2] = -0.5
	P[6] = -0.5  
	return net, E, X, P 
end


function sourcesinklocs(ns::Int64, nd::Int64, n::Int64)
	sources = Array{Int64, 1}([])
	sinks = Array{Int64, 1}([])
	for i=1:ns
		while true
			z = rand(1:n)
			if z ∉ sources
				push!(sources, z)
				break 
			end
		end
	end
	for i=1:nd
		while true
			z = rand(1:n)
			if (z ∉ sinks) & (z ∉ sources) 
				push!(sinks, z)
				break 
			end 
		end
	end
	return sources, sinks 
end


function sourcesinklocsmod(net::Graph, ns::Int64, nd::Int64, n::Int64)
	sources = Array{Int64, 1}([])
	sinks = Array{Int64, 1}([])
	if ns == 1 
		while true
			z = rand(1:n)
			if length(net.nodelist[z].outedges) + length(net.nodelist[z].inedges) > 3 
				push!(sources, z)
				break 
			end
		end 
	end
	if nd == 1
		while true 
			z = rand(1:n)  
			if length(net.nodelist[z].outedges) + length(net.nodelist[z].inedges) > 3
				push!(sinks, z)
				break
			end
		end
	end
	if ns > 1	
		for i=1:ns
			while true
				z = rand(1:n)
				if (z ∉ sinks) & (z ∉ sources) 
					push!(sources, z)
					break 
				end
			end
		end
	end 
	if nd > 1
		for i=1:nd
			while true
				z = rand(1:n)
				if (z ∉ sinks) & (z ∉ sources) 
					push!(sinks, z)
					break 
				end 
			end
		end
	end
	return sources, sinks 
end 


function sourcesinkwellmixed(ns::Int64, nd::Int64, n::Int64)
	sources = Array{Int64, 1}([])
	sinks = Array{Int64, 1}([])
	for i=1:2:n 
		push!(sources, i)
	end
	for i=1:n 
		if i ∉ sources 
			push!(sinks, i)
		end
	end 
	return sources, sinks 
end


function sourcesinkhalfhalf(ns::Int64, nd::Int64, n::Int64)
	sources = Array{Int64, 1}([])
	sinks = Array{Int64, 1}([])
	for i=1:div(n,2)
		push!(sources, i)
	end
	for i=1:n 
		if i ∉ sources 
			push!(sinks, i)
		end
	end
	return sources, sinks 
end


function sourcessinkclockface(n::Int64, pos::Int64)
	sources = Array{Int64, 1}([])
	sinks = Array{Int64, 1}([])
	push!(sources,1)
	push!(sources, pos)
	for i=2:n
		if i ≠ pos 
			push!(sinks, i)
		end
	end
	return sources, sinks 
end


function sourcesinkvector(sources::Vector{Int64}, sinks::Vector{Int64}, n::Int64, Ptot::Float64)
	P = zeros(Float64, n)
	sourcenum = length(sources)
	sinknum = length(sinks)
	for i ∈ sources
		P[i] = Ptot/Float64(sourcenum)
	end
	for i ∈ sinks
		P[i] = -Ptot/Float64(sinknum) 
	end
	return P
end


function makenoise(sourcenum::Int64, sinknum::Int64) 
	d = Normal(0.0, 0.4) 
	sourcenoise = rand(d, sourcenum)
	sinknoise = rand(d, sinknum) 
	return sourcenoise, sinknoise 
end 


function makegammasourcesinkvector(sources::Vector{Int64}, sinks::Vector{Int64}, n::Int64, Ptot::Float64)
	P = P = zeros(Float64, n)
	sourcenum = length(sources)
	sinknum = length(sinks)
	d = Gamma(2.0, 0.5) 
	runningsum = 0.0 
	for i ∈ sources 
		P[i] = rand(d) 
		runningsum += P[i] 
	end
	alpha = Ptot/runningsum 
	for i ∈ sources
		P[i] *= alpha
	end
	runningsum = 0.0 
	for i ∈ sinks 
		P[i] = -rand(d) 
		runningsum += abs(P[i])
	end
	alpha = Ptot/runningsum
	for i ∈ sinks 
		P[i] *= alpha
	end
	return P
end	


function noisysourcesinkvector(sources::Vector{Int64}, sinks::Vector{Int64}, n::Int64, Ptot::Float64)
	P = zeros(Float64, n)
	sourcenum = length(sources)
	sinknum = length(sinks)
	sourcenoise, sinknoise = makenoise(sourcenum, sinknum)  
	runningsum = 0.0 
	for (i, val) ∈ enumerate(sources)  
		P[val] = (Ptot/Float64(sourcenum))*(1.0 + sourcenoise[i])
		runningsum += P[val]
	end
	delta = Ptot - runningsum 
	for val ∈ sources 
		P[val] += delta/sourcenum 
	end 
	runningsum = 0.0  
	for (i, val) ∈ enumerate(sinks)  
		P[val] = -(Ptot/Float64(sinknum))*(1.0 + sinknoise[i])
		runningsum += P[val]
	end
	delta = -Ptot - runningsum 
	for val ∈ sinks 
		P[val] += delta/sinknum  
	end 
	return P 
end


















