using .NetworkResilience
using Statistics 

function main(t)

	ind = parse(Int64, ARGS[1])

	configs = [] 
	n = 40
	for j=1:n
		for i=1:n-j 
			ne = n - i - j 
			push!(configs, [i j ne])
		end
	end
	break1 = length(configs) 
	# println(break1)

	# n = 100  
	# for j=1:n
	# 	for i=1:n-j 
	# 		ne = n - i - j 
	# 		push!(configs, [i j ne])
	# 	end
	# end
	# # break2 = length(configs) 
	# println(break1)
	# # println(break2-break1)
	# # println(break2)

	# # # n = 64
	# # # for j=1:n
	# # # 	for i=1:n-j 
	# # # 		ne = n - i - j 
	# # # 		push!(configs, [i j ne])
	# # # 	end
	# # # end
	# # # break3 = length(configs)

	#Make the file name 
	directory = "redoneq0/"
	indexdict = Dict("40" => ind, "100" => ind-break1)
	subdirectorydict = Dict("40" => "40tri/", "100" => "100tri/")
	ntot = sum(configs[ind])
	thisid = indexdict[string(ntot)]
	thissubdir = subdirectorydict[string(ntot)]
	fname = "data/" * directory * thissubdir * string(thisid) * ".txt"

	#Compute the cascade 
	ns = configs[ind][1]
	nd = configs[ind][2]
	ne = configs[ind][3]
	ensemblesize = 200 
	# subensemblesize = 40 
	tol = 5e-4 
	n = ns + nd + ne 
	q = 0.0
	k = 4 
	Ptot = 1.0 
	# ac, fmax, iter = cascadebisectionrandomtrigger(ensemblesize, subensemblesize, tol, n, ns, nd, q, k, Ptot) 
	ac, fmax, iter = cascadebisection2(ensemblesize, tol, n, ns, nd, q, k, Ptot)

	open(fname, "w") do f 
		write(f, "$ns $nd $ne $fmax \n")
		write(f, "$ac $fmax $iter $iter \n")
	end

	# ensemblesize = 1 
	# topensemblesize = 2000
	# q = 0.1
	# alphares = 400
	# alphamin = 0.2
	# alphamax = 2.0
	# k = 4 
	# Ptot = 1.0
	# Xval = 1.0
	# es = zeros(Float64, topensemblesize, alphares)
	# ts = zeros(Float64, topensemblesize, alphares)
	# as = zeros(Float64, alphares)
	# fmaxs = zeros(Float64, topensemblesize)
	# for i=1:topensemblesize

	# 	a, e, std, t, tstd, mate, matt, fmax = smcascadefixedtop(ensemblesize, ntot, ns, nd, q, alphares, alphamin, alphamax, k, Ptot, Xval)
	# 	es[i,:] = e 
	# 	ts[i,:] = t 
	# 	as[:] = a
	# 	fmaxs[i] = fmax

	# end 

	# #Compute mean over the topolgy ensemble 
	# eav = mean(es, dims=1)[1,:]
	# tav = mean(ts, dims=1)[1,:] 
	# estd = std(es, dims=1)[1,:]
	# meanfmax = mean(fmaxs) 

	#Write stuff to the file
	# open(fname, "w") do f 
	# 	write(f, "$ns $nd $ne $meanfmax \n")
	# 	for i=1:alphares
	# 		beta = eav[i]
	# 		gamma = estd[i]
	# 		kappa = as[i]
	# 		epsilon = tav[i]
	# 		write(f, "$kappa $beta $gamma $epsilon \n")
	# 	end
	# end

end

main(ARGS)
















