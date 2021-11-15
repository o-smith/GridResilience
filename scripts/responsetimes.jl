using .NetworkResilience
using Statistics

ensemblesize = 200
n = 50 
ns = 10 
nd = 40 
q = 0.1 
k = 4 
alphares = 50 
alphamin = 0.5 
alphamax = 2.5 
I = 1.0 
kappa = 5.0 

println("Computing cascade durations for a range of gamma values...")
Ds = [0.3:0.1:2.0;]
for D in Ds 
	D_alias = Int(10*D)
	fname = "data/durations/q1_10_40_0_t_gamma$D_alias.txt" 

	a, s, t, fmax, crits = swingcascadeandtvalpha(ensemblesize, n, ns, nd, q, k, alphares, alphamin, alphamax, I, D, kappa)
	open(fname, "w") do f 
		for i=1:length(a)
			beta = a[i]
			gamma = t[i]
			write(f, "$beta $gamma \n")
		end
	end
end 


D = 1.0 
q = 0.0
qs = [0.0:0.1:1.0;]
println("Computing cascade durations for a range of q values...")
for q in qs 
	q_alias = Int(10*q)
	fname = "data/durations/q1_10_40_0_t_q$q_alias.txt" 

	a, s, t, fmax, crits = swingcascadeandtvalpha(ensemblesize, n, ns, nd, q, k, alphares, alphamin, alphamax, I, D, kappa)
	open(fname, "w") do f 
		for i=1:length(a)
			beta = a[i]
			gamma = t[i]
			write(f, "$beta $gamma \n")
		end
	end
end 


q = 0.1
kappa = 0.3
kappas = [1.0:1.0:11.0;]
println("Computing cascade durations for a range of kappa values...")
for kappa in kappas 
	kappa_alias = Int(kappa)
	fname = "data/durations/q1_10_40_0_t_kappa$kappa_alias.txt" 

	a, s, t, fmax, crits = swingcascadeandtvalpha(ensemblesize, n, ns, nd, q, k, alphares, alphamin, alphamax, I, D, kappa)
	open(fname, "w") do f 
		for i=1:length(a)
			beta = a[i]
			gamma = t[i]
			write(f, "$beta $gamma \n")
		end
	end
end 


kappa = 5.0 
println("Computing cascade durations for a range of network sizes n...")
ns = [45:5:100;] 
for n in ns
	ns = Int(n/5) 
	nd = Int(4*n/5)
	fname = "data/durations/q1_t_$n.txt"

	a, s, t, fmax, crits = swingcascadeandtvalpha(ensemblesize, n, ns, nd, q, k, alphares, alphamin, alphamax, I, D, kappa)
	open(fname, "w") do f 
		for i=1:length(a)
			beta = a[i]
			gamma = t[i]
			write(f, "$beta $gamma \n")
		end
	end

end
