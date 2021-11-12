using .NetworkResilience
using PyPlot
using Statistics
using Distributed

#Define colours 
martaRed = "#c24c51" 
martaGreen = "#54a666"
martaBlue = "#4c70b0"
martaPurple = "#7f70b0"
martaGold = "#ccb873"
rc("font", family="serif", size=20)


net, E = austriannetwork()
println(numnodes(net))
println(numedges(net))
println(mean(getdegrees(E)))


# net, E = austriannetwork()
# A = adjacencymatrix(E)
# n = numnodes(net)
# m = numedges(net)
# println(n, "    ", m)
# components, table = connectedcomponents(A)
# println(components)

# chunk = table[1]
# println(chunk)
# nodeset = collect(1:n)
# nodesubset = zeros(Int64, length(chunk))
# for j=1:length(chunk)
# 	nodesubset[j] = nodeset[chunk[j]]
# end
# println(nodesubset)

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
# println(edgecardinality)
# nodecardinality = length(nodeset)

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

# imshow(E1)
# show() 
# A = adjacencymatrix(E1)
# components, table = connectedcomponents(A)
# println(components)

# # println(n, "    ", m)
# # imshow(A, interpolation="none", cmap="seismic")
# # show() 



ntot, ns, nd = 67, 2, 65
ensemblesize = 1
topensemblesize = 100
q = 0.1
alphares = 150
alphamin = 0.01
alphamax = 2.0
k = 4 
Ptot = 1.0
Xval = 1.0
es = zeros(Float64, topensemblesize, alphares)
ts = zeros(Float64, topensemblesize, alphares)
as = zeros(Float64, alphares)
fmaxs = zeros(Float64, topensemblesize)
for i=1:topensemblesize
	println(i)
	a, e, std, t, tstd, mate, matt, fmax = smcascadefixedtop(ensemblesize, ntot, ns, nd, q, alphares, alphamin, alphamax, k, Ptot, Xval)
	es[i,:] = e 
	ts[i,:] = t 
	as[:] = a
	fmaxs[i] = fmax
end 
eav = mean(es, dims=1)[1,:]
tav = mean(ts, dims=1)[1,:] 
estd = std(es, dims=1)[1,:]
meanfmax = mean(fmaxs) 
aftertransition = findall(x->x>0.5, eav)
# ac = as[aftertransition[1]]	
plot(as, eav, "-", lw=3.0, color=martaBlue)
xlabel("\$\\alpha\$", fontsize=24)
ylabel("\$\\mathcal{S}\$", rotation=0, labelpad=0, fontsize=24)
# xlim([0,2.5])
# xticks([0,1])
ylim([0,1])
yticks([0,1])
# axvline(meanfmax, linestyle="--", linewidth=3.0, color="k")
# axvline(ac, linestyle="--", linewidth=5.0, color="k")
tight_layout()
println(meanfmax)
# println(ac)
# println(ac/meanfmax)
show() 

# fname = "data/10_40_0_0p1_200profile.txt"
# open(fname, "w") do f 
# 	write(f, "$meanfmax $ac \n")
# 	for i=1:100
# 		beta = as[i] 
# 		gamma = eav[i] 
# 		write(f, "$beta $gamma \n")
# 	end 
# end

# n = 200 
# ns = 100 
# nd = 100 
# q = 0.1 
# k = 4 
# Ptot = 1.0 
# net, E = smallworldnetwork(n, k, q) 
# sources, sinks = sourcesinklocs(ns, nd, n)
# m = numedges(net)
# X = zeros(Float64, m) .+ 1.0
# P = sourcesinkvector(sources, sinks, n, Ptot)
# L = weightedlaplacianmatrix(E, X) 
# theta = nodevoltages(L, P)
# F = edgecurrents(theta, E, X)
# probs = triggerhist(F) 
# hist(probs)
# show() 

# ac, fmax, iter = cascadebisectionrandomtrigger(100, 20, 5e-4, 40, 10, 30, 0.1, 4, 1.0)
# println(ac)

# ac, fmax, iter = cascadebisection(100, 5e-4, 14, 2, 12, 0.1, 4, 1.0)
# println(ac)


# n = 32
# as = []
# ks = [] 
# for k=2:2:(n+1)
# 	println(k)
# 	a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	aftertransition = findall(x->x>0.5, e)
# 	ac = a[aftertransition[1]] - fmax
# 	push!(as, ac)
# 	push!(ks, k)
# end 
# as1 = as
# ks1 = ks


# plot(ks, as, ".-", color=martaRed, lw=3.0, label="\$n=32\$")


# n = 24
# as = []
# ks = [] 
# for k=2:2:(n+1)
# 	println(k)
# 	a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	aftertransition = findall(x->x>0.5, e)
# 	ac = a[aftertransition[1]] - fmax
# 	push!(as, ac)
# 	push!(ks, k)
# end 
# plot(ks, as, ".-", color=martaBlue, lw=3.0, label="\$n=24\$")

# n = 16
# as = []
# ks = [] 
# for k=2:2:(n+1)
# 	println(k)
# 	a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	aftertransition = findall(x->x>0.5, e)
# 	ac = a[aftertransition[1]] - fmax
# 	push!(as, ac)
# 	push!(ks, k)
# end 
# plot(ks, as, ".-", color=martaGreen, lw=3.0, label="\$n=16\$")
# plot(ks1, 1.0./(ks1.*(ks1.-1.0)), "k--", lw=3.0, label="\$\\frac{P}{k(k-1)}\$")
# ylabel("\$\\widetilde{\\alpha}\$", rotation=0, labelpad=20)
# xlabel("\$k\$")
# legend()
# tight_layout() 
# show()


# #tentative simplex code ####################
# n = 32
# mat = zeros(Float64, n, n)
# for j=1:n
# 	for i=1:n-j
# 		println(i," ", j)
# 		a, e, std, t, tstd, p0, p1, p2, mato = smcascadefixedtop(100, n, i, j, 0.0, 100, 0.01, 10.0, 4, 1.0, 1.0)
# 		aftertransition = findall(x->x>0.5, e)
# 		ac = a[aftertransition[1]]
# 		# println(ac)
# 		mat[i,j] = ac 
# 	end
# 	println("poop")
# 	println(j)
# end

# for i=1:n 
# 	for j=1:n 
# 		if mat[i,j] == 0.0
# 			mat[i,j] = NaN
# 		end
# 	end
# end
# imshow(mat, cmap="viridis", interpolation="none")
# colorbar()
# show()
# ############################
# break1 = len(configs) 

# a, e, std, t, tstd, p0, p1, p2, mat = smcascadefixedtop(100, n, ns, nd, 0.0, 100, 0.01, 2.0, 4, 1.0, 1.0)
# aftertransition = findall(x->x>0.5, e)
# ac = a[aftertransition[1]]


# plot(ns, ac, "bo")
# plot(a, e, color="k", ls="-", lw=3.0)
# axvline(a[10], color="r")
# xlabel("Alpha")
# ylabel("Remaining edges")
# show() 


# fill_between(a, t.+tstd, t.-tstd, color="k", alpha=0.2)
# plot(a, t, color="k", ls="-", lw=3.0)
# xlabel("α")
# ylabel("Time taken")
# show() 

# for j=1:size(mat,1)
# 	plot(a, mat[j,:], "b.")
# end
# ylim([0,1])
# show()

# hist(p0, 50, alpha=0.75)
# show()

# hist(p1, 50, alpha=0.75)
# show()

# hist(p2, 50, alpha=0.75)
# show()

# edgemat = [] 
# for i=1:2000
# 	println(i)
# 	ralpha, rT, redges, rsources, rweightedsources, rsinks, rpower = smallworldalphas(64,32,32,0.0, 0.02, 1.0, 10.0, 4, 20000.0, 0.5)
# 	if i==1
# 		global edgemat = redges
# 	else
# 		vcat(edgemat, redges)
# 	end
# 	plot(ralpha, redges)
# end 
# show() 
# println(length(edgemat))
# ed = mean(edgemat, dims=2)
# plot(ralpha, ed)
# show() 





