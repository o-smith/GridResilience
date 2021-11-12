using .NetworkResilience
using Statistics 

function main(t)

	ind = parse(Int64, ARGS[1])

	res = 40
	configs = [] 
	qs = collect(range(0, length=res, stop=1.0)) 
	ns = collect(range(42, length=res, stop=276)) 
	for i=1:res
		for j=1:res
			push!(configs, [i j qs[i] ns[j]])
		end
	end

	#Make the file name 
	directory = "data/qntrans_50_50/"
	fname = directory * string(ind) * ".txt" 

	#Set some parameters 
	Ptot = 1.0
	k = 4
	ensemblesize = 2500 
	tol = 5e-4 
	n = Int64(configs[ind][4]) 
	q = configs[ind][3] 

	#Compute cascade 
	ac_corner, fmax_corner, iter, degss = cascadebisection2(ensemblesize, tol, n, 1, n-1, q, k, Ptot)
	ac_middle, fmax_middle, iter, degss = cascadebisection2(ensemblesize, tol, n, div(n,6), div(n,6), q, k, Ptot)

	#Save data 
	i = Int64(configs[ind][1])
	j = Int64(configs[ind][2]) 
	open(fname, "w") do f 
		write(f, "$i $j $q $n $ac_corner $ac_middle \n")
	end

end

main(ARGS)
















