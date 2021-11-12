using .NetworkResilience
using Statistics 

function main(t)

	ind = parse(Int64, ARGS[1])

	configs = [] 
	for i=10:2:120
		ns = div(i,2)
		push!(configs, [ns ns 0])
	end
	for i=120:10:400
		ns = div(i,2)
		push!(configs, [ns ns 0])
	end 

	#Make the file name 
	directory = "tri_sm_1/nconthalfhalf/"
	fname = "data/" * directory * string(ind) * ".txt"

	#Compute the cascade 
	ns = configs[ind][1]
	nd = configs[ind][2]
	ne = configs[ind][3] 
	ntot = sum(configs[ind])
	ensemblesize = 20 
	topensemblesize = 20
	q = 0.1
	alphares = 300
	alphamin = 0.01
	alphamax = 0.5
	k = 4 
	Ptot = 1.0
	Xval = 1.0
	es = zeros(Float64, topensemblesize, alphares)
	ts = zeros(Float64, topensemblesize, alphares)
	as = zeros(Float64, alphares)
	fmaxs = zeros(Float64, topensemblesize)
	for i=1:topensemblesize

		a, e, std, t, tstd, mate, matt, fmax = smcascadefixedtop(ensemblesize, ntot, ns, nd, q, alphares, alphamin, alphamax, k, Ptot, Xval)
		es[i,:] = e 
		ts[i,:] = t 
		as[:] = a
		fmaxs[i] = fmax

	end 

	#Compute mean over the topolgy ensemble 
	eav = mean(es, dims=1)[1,:]
	tav = mean(ts, dims=1)[1,:] 
	estd = std(es, dims=1)[1,:]
	meanfmax = mean(fmaxs) 

	#Write stuff to the file
	open(fname, "w") do f 
		write(f, "$ns $nd $ne $meanfmax \n")
		for i=1:alphares
			beta = eav[i]
			gamma = estd[i]
			kappa = as[i]
			epsilon = tav[i]
			write(f, "$kappa $beta $gamma $epsilon \n")
		end
	end

end

main(ARGS)













