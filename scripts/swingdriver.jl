using .NetworkResilience
using Statistics 
using PyPlot
# using DifferentialEquations 
# using LinearAlgebra 

ensemblesize = 200
tol = 5e-4 
n = 60
ns = 30
nd = 30
q = 1.0
I = 1.0 
D = 1.0 
k = 4 
Ptot = 1.0 
κ = 5.0 
alphas, s, ov, desy, unb, fmaxes, crit = swingcascadebreakdownvsalpha(ensemblesize, n, ns, nd, q, k, 50, 0.1, 2.5, I, D, κ)  

plot(alphas, s, color=:blue)
plot(alphas, ov, color=:red)
plot(alphas, desy, color=:orange)
plot(alphas, unb, color=:black) 
show() 
# println(α) 

# alphas, S, fmax, critvalues = swingcascadevsalpha(ensemblesize, n, ns, nd, q, 4, 50, 0.1, 2.5)

# aftertransition = findall(x->x>0.5, s)
# ac = alphas[aftertransition[1]]	
# println(ac)

# # fname = "data/new/q1_30_30_0_profile.txt"
# # open(fname, "w") do f 
# # 	for i=1:length(alphas) 
# # 		a = alphas[i]
# # 		b = s[i]
# # 		write(f, "$a $b \n")
# # 	end 
# # end 

# fname = "data/new/q0_30_30_0.txt"
# open(fname, "w") do f 
# 	for i=1:length(alphas) 
# 		a = alphas[i]
# 		b = s[i]
# 		c = ov[i]
# 		d = desy[i] 
# 		e = unb[i]
# 		write(f, "$a $b $c $d $e \n")
# 	end 
# end 

# plot(alphas, s)
# savefig("prof.pdf")
# x = mean(filter(!isnan, rhos))
# println(x)
# # # hist(critvalues)
# hist(rhos)
# # axvline(ac, linestyle="--", linewidth=5.0, color="k")
# show() 



