using .NetworkResilience
using Statistics 


configs = [] 
n = 30
for j=1:n
	for i=1:n-j 
		ne = n - i - j 
		push!(configs, [i j ne])
	end
end
#435 points 

ensemblesize = 50
q = 0.1 
alphares = 500 
alphamin = 0.01 
alphamax = 0.5 
k = 4 
Ptot = 1.0 
Xval = 1.0 

#Create a network topology
net, E0 = smallworldnetwork(n, k, q)

for z=200:350
	println(z) 

	ns = configs[z][1]
	nd = configs[z][2]
	ne = configs[z][3]  
	a, e, std, t, tstd, mate, matt, fmax = singlecascade(net, E0, ensemblesize, n, ns, nd, q, alphares, alphamin, alphamax, k, Ptot, Xval)

	#Make the file name 
	directory = "tri_sm_1_single2/"
	ntot = sum(configs[z])
	fname = "data/" * directory * string(z) * ".txt" 

	#Write stuff to the file
	open(fname, "w") do f 
		write(f, "$ns $nd $ne $fmax \n")
		for i=1:alphares
			beta = e[i]
			gamma = std[i]
			kappa = a[i]
			epsilon = t[i] 
			write(f, "$kappa $beta $gamma $epsilon \n")
		end
	end

end 

