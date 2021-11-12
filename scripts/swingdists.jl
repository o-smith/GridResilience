using .NetworkResilience 
using PyPlot 
using Statistics 



n = 100
ns = 1
nd = 99 
q = 1.0
κ = 5.0
I = 1.0
D = 1.0 
tol = 5e-4
k = 4
ensemblesize = 200
alphamin=0.1
alphamax=2.5
alphares=50
acs = [] 

# ac, fmax, iter = cascadebisection2(ensemblesize, tol, n, ns, nd, q, k, 1.0)
# println(ac) 

alp, s, t, fmax, crits = swingcascadeandtvalpha(ensemblesize, n, ns, nd, q, k,alphares,alphamin,alphamax,
								I, D, κ)
for i=1:length(alp)
	println(alp[i], ": ", s[i], ": ", t[i])
end

plot(alp,t)
savefig("im2.pdf")

# fname = "data/new/q1_15_45_0_time.txt"
# open(fname, "w") do f 
# 	for i=1:length(alp)
# 		println(i)
# 		x = alp[i]
# 		y = s[i]
# 		z = t[i] 
# 		g = crits[i]
# 		write(f, "$x $y $z $g \n") 
# 	end
# end


# n = 50
# ns = 10
# nd = 40
# q = 0.1
# κ = 5.0
# I = 1.0
# D = 1.0 
# tol = 5e-4
# k = 4
# ensemblesize = 500
# alpha = 2.0
# acs = []

# s, t = svt(ensemblesize, n, ns, nd, q, k,alpha, I, D, κ)


# fname = "data/swingtimeprofs/q1_10_40_0_svt_alpha2.txt"
# open(fname, "w") do f 
# 	for i=1:length(s)
# 		println(i)
# 		x = s[i]
# 		y = t[i]
# 		write(f, "$x $y  \n") 
# 	end
# end



# for i=1:ensemblesize
# 	println(i) 
# 	ac, fmax, iter = swingcascadebisection(1, tol, n, ns, nd, q, k, κ, I, D)
# 	push!(acs, ac) 
# end


# hist(acs)
# savefig("im4.pdf")

# fname = "data/new/hist_q1_30_30_0.txt"
# open(fname, "w") do f 
# 	for i=1:length(acs)
# 		x = acs[i]
# 		write(f, "$x \n") 
# 	end
# end