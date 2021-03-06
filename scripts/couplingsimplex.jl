using .NetworkResilience
using Statistics 

function main(t)

	ind = parse(Int64, ARGS[1])

	configs = [] 
	n = 67
	for j=1:n
		for i=1:n-j 
			ne = n - i - j 
			push!(configs, [i j ne])
		end
	end
	# break1 = length(configs) 

	# #Make the file name 
	# directory = "swing_coupling_q0/"
	# indexdict = Dict("50" => ind, "100" => ind-break1)
	# subdirectorydict = Dict("50" => "50tri/", "100" => "100tri/")
	# ntot = sum(configs[ind])
	# thisid = indexdict[string(ntot)]
	# thissubdir = subdirectorydict[string(ntot)]
	# fname = "data/" * directory * thissubdir * string(thisid) * ".txt"

	#Compute the cascade 
	ns = configs[ind][1]
	nd = configs[ind][2]
	ne = configs[ind][3]
	ensemblesize = 3
	tol = 5e-4 
	n = ns + nd + ne 
	q = 0.0
	I = 1.0 
	D = 1.0 
	k = 4
	Ptot = 1.0 
	κ = 5.5
	criticalcouplings = criticalcouplingbisection(ensemblesize, tol, n, ns, nd, q, k, κ, I, D)
	println(criticalcouplings)
	# open(fname, "w") do f 
	# 	for z=1:ensemblesize
	# 		x = criticalcouplings[z]
	# 		write(f, "$x \n")
	# 	end
	# end

end

main(ARGS)
















