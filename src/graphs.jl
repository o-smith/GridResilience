mutable struct Node
	nodeindex::Int64
	outedges::Array{Int64,1}
	inedges::Array{Int64,1} 
end

mutable struct Edge
	edgeindex::Int64
	origin::Int64
	target::Int64
	active::Bool 
end

mutable struct Graph
	nodelist::Array{Node,1}
	edgelist::Array{Edge,1} 
	edgecounter::Int64
	nodecounter::Int64 
end 

function addnode!(graph::Graph)
	i = graph.nodecounter
	v = Node(i, Int64[], Int64[]) 
	push!(graph.nodelist, v) 
	graph.nodecounter += 1  
end 

function addedge!(graph::Graph, start::Int64, fin::Int64)
	i = graph.edgecounter
	e = Edge(i, start, fin, true) 
	push!(graph.edgelist, e) 
	push!(graph.nodelist[start].outedges, i) 
	push!(graph.nodelist[fin].inedges, i)
	graph.edgecounter += 1 
end 

function numedges(graph::Graph)
	return length(graph.edgelist)
end

function numnodes(graph::Graph)
	return length(graph.nodelist)
end 

mutable struct NetStore
	adjmat::Matrix{Float64}
	pvec::Vector{Float64} 
	kappa::Float64
end 

mutable struct NetStore2
	adjmat::Matrix{Float64}
	pvec::Vector{Float64} 
	incidence::Matrix{Float64}
	kappa::Float64
	synctol::Float64
	alpha::Float64
	maxflow::Float64
	sync::Bool
	trip::Bool 
end 







