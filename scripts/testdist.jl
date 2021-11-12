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

n = 20 
k = 4 
q = 0.0 
ns = 10 
nd = 10 
Ptot = 1.0 
net, E = smallworldnetwork(n, k, q) 
sources, sinks = sourcesinklocs(ns, nd, n)
P = makegammasourcesinkvector(sources, sinks, n, Ptot)
println(P)
println(sum(P))

# ntot, ns, nd = 100, 1, 99
# ensemblesize = 100 #ntot*2
# topensemblesize = 100
# q = 0.1
# alpha = 1.8
# k = 4 
# Ptot = 1.0
# Xval = 1.0
# rloads = zeros(Float64, ensemblesize*topensemblesize) 
# redges = zeros(Float64, ensemblesize*topensemblesize) 
# rpower = zeros(Float64, ensemblesize*topensemblesize)  

# for i=1:topensemblesize
# 	println(i)
# 	relativeload, survivingedges, survivingpower = smcascadedist(ensemblesize, ntot, ns, nd, q, alpha, k, Ptot, Xval)
# 	rloads[(i-1)*ensemblesize+1:i*ensemblesize] = relativeload
# 	redges[(i-1)*ensemblesize+1:i*ensemblesize] = survivingedges
# 	rpower[(i-1)*ensemblesize+1:i*ensemblesize] = survivingpower
# end 
# # hist2d(rloads, redges) #, "o", color=martaRed) 
# # show() 

# #Write stuff to the file
# fname = "testdata.txt"
# open(fname, "w") do f 
# 	# write(f, "$ns $nd $ne \n")
# 	for i=1:length(rloads)
# 		beta = rloads[i]
# 		gamma = redges[i]
# 		kappa = rpower[i]
# 		write(f, "$beta $gamma $kappa \n")
# 	end
# end
# ensemblesize = 100
# n = 64
# ns = 10
# nd = 10
# q = 1.0
# alphares = 500
# alphamin = 0.01
# alphamax = 0.4
# k = 4 
# Ptot = 1.0
# Xval = 1.0

# topensemblesize = 1
# es = zeros(Float64, topensemblesize, alphares)
# as = zeros(Float64, alphares)
# acs = zeros(Float64, topensemblesize)
# fmaxs = zeros(Float64, topensemblesize)
# for i=1:topensemblesize
# 	println("topen = ", i)

# 	a, e, std, t, tstd, mate, matt, fmax = smcascadefixedtop(ensemblesize, n, ns, nd, q, alphares, alphamin, alphamax, k, Ptot, Xval)
# 	es[i,:] = e 
# 	as[:] = a
# 	fmaxs[i] = fmax

# 	aftertransition = findall(x->x>0.5, e)
# 	ac = a[aftertransition[1]]
# 	println(ac)
# 	acs[i] = ac

# 	for j=1:size(mate,1)
# 		plot(a, mate[j,:], ".", color=martaBlue)
# 	end

# end 


# println("Fmax av = ", mean(fmaxs), " for q = ", q, " ns= ", ns)
# eav = mean(es, dims=1)[1,:]
# estd = std(es, dims=1)[1,:]

# plot(as, eav, "o", color=martaRed)
# axvline(mean(fmaxs), color="k", alpha=0.4, lw=2.0)
# fill_between(as, eav.+estd, eav.-estd, color=martaRed, alpha=0.2)

# xlabel("\$\\alpha\$")
# ylabel("\$ {\\mathcal{S}}\$", rotation=0, labelpad=20)
# ylim([0,1.2])
# tight_layout() 
# show() 

# hist(acs, color=martaRed)
# # ylabel("\$P(\\alpha_c)\$", rotation=0, labelpad=20)
# xlabel("\$\\alpha_c\$")
# tight_layout() 
# show() 






# a, e, std, t, tstd, mate, matt, fmax = smcascadefixedtop(ensemblesize, n, ns, nd, q, alphares, alphamin, alphamax, k, Ptot, Xval)



# plot(a, e, "o", color=martaRed)
# fill_between(a, e.+std, e.-std, color=martaRed, alpha=0.2)

# for j=1:size(mate,1)
# 	plot(a, mate[j,:], ".", color=martaBlue)
# end

# xlabel("\$\\alpha\$")
# ylabel("\$ {\\mathcal{S}}\$", rotation=0, labelpad=20)
# ylim([0,1.2])
# tight_layout() 
# show() 































