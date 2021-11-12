using .NetworkResilience
using PyPlot
using Statistics

#Define colours 
martaRed = "#c24c51" 
martaGreen = "#54a666"
martaBlue = "#4c70b0"
martaPurple = "#7f70b0"
martaGold = "#ccb873"
rc("font", family="serif", size=20)

n = 40
as = []
ks = [] 
# for k=4:2:5
# 	println(k)
# 	# a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	# aftertransition = findall(x->x>0.5, e)
# 	# ac = a[aftertransition[1]]
# 	# println(fmax, 1.0/k)
# 	ac, fmax, iter = cascadebisection2(500, 5e-5, n, 1, 39, 1., k, 1.0)
# 	push!(as, ac)
# 	push!(ks, k)
# end 
# plot(ks, as, ".", color=martaBlue, ms=12.0, label="\$n=20\$")

for k=1:20
	println(k)
	# a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
	# aftertransition = findall(x->x>0.5, e)
	# ac = a[aftertransition[1]]
	# println(fmax, 1.0/k)
	ac, fmax, iter = cascadebisection(200, 5e-5, 40, k, k, 0.2, 4, 1.0)
	push!(as, ac)
	# push!(ks, k)
end 
plot(as, ".", color=martaBlue, ms=12.0, label="\$n=20\$")

# n = 30
# as = []
# ks = [] 
# for k=2:2:15
# 	println(k)
# 	# a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	# aftertransition = findall(x->x>0.5, e)
# 	# ac = a[aftertransition[1]]
# 	# println(fmax, 1.0/k)
# 	ac, fmax, iter = cascadebisection(1, 5e-5, n, 1, n-1, 0.0, k, 1.0)
# 	push!(as, ac)
# 	push!(ks, k)
# end 
# plot(ks, as, ".", color=martaGreen, ms=12.0, label="\$n=30\$")

# n = 40
# as = []
# ks = [] 
# for k=2:2:15
# 	println(k)
# 	# a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	# aftertransition = findall(x->x>0.5, e)
# 	# ac = a[aftertransition[1]]
# 	# println(fmax, 1.0/k)
# 	ac, fmax, iter = cascadebisection(1, 5e-5, n, 1, n-1, 0.0, k, 1.0)
# 	push!(as, ac)
# 	push!(ks, k)
# end 
# plot(ks, as, ".", color=martaRed, ms=12.0, label="\$n=40\$")
# # fname = "data/kcont/n64_q0.txt"
# open(fname, "w") do f 
# 	for i=1:length(as)
# 		kdat = ks[i]
# 		rho = as[i]
# 		write(f, "$kdat $rho \n")
# 	end
# end 

# x = 1.01
# curve = (ks)./(ks.-1.0) 
# # plot(ks, as, ".-", color=martaGreen, lw=3.0, label="\$n=16\$")
# plot(ks, curve, "--", color="k", lw=3.0, label="\$K/(K-1)\$")
xlabel("\$K\$")
ylabel("\$\\rho\$", rotation=0,labelpad=20)
# legend()
# plot(ks, ks./(ks.-1), ".-", color="b", lw=3.0, label="\$n=16\$")
# plot(ks, 2.0.*ks./(3.0.*(ks.-1)), ".-", color="m", lw=3.0, label="\$n=16\$")


# n = 24
# as = []
# ks = [] 
# kmns = [] 
# for k=4:2:40
# 	# a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	# aftertransition = findall(x->x>0.5, e)
# 	# ac = a[aftertransition[1]]
# 	# println(fmax, 1.0/k)
# 	# println(k)
# 	ac, fmax, iter, degmeans = cascadebisection(1, 5e-5, n, 1, n-1, 0.0, k, 1.0)
# 	push!(as, ac)
# 	push!(ks, k)
# 	push!(kmns, mean(degmeans)) 
# end 
# fname = "data/kcont/n24_q0.txt"
# open(fname, "w") do f 
# 	for i=1:length(as)
# 		kdat = ks[i]
# 		rho = as[i]
# 		write(f, "$kdat $rho \n")
# 	end
# end 

# n = 32
# as = []
# ks = [] 
# kmns = [] 
# for k=4:2:40
# 	# a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	# aftertransition = findall(x->x>0.5, e)
# 	# ac = a[aftertransition[1]]
# 	# println(fmax, 1.0/k)
# 	# println(k)
# 	ac, fmax, iter, degmeans = cascadebisection(1, 5e-5, n, 1, n-1, 0.0, k, 1.0)
# 	push!(as, ac)
# 	push!(ks, k)
# 	push!(kmns, mean(degmeans)) 
# end 
# fname = "data/kcont/n32_q0.txt"
# open(fname, "w") do f 
# 	for i=1:length(as)
# 		kdat = ks[i]
# 		rho = as[i]
# 		write(f, "$kdat $rho \n")
# 	end
# end 

# fname = "data/kcont/n70_q5.txt"
# open(fname, "w") do f 
# 	for i=1:length(as)
# 		beta = ks[i]
# 		gamma = as[i]
# 		kappa = kmns[i] 
# 		write(f, "$beta $gamma $kappa \n")
# 	end
# end

# plot(ks, as, ".-", color=martaBlue, lw=3.0, label="\$n=32\$")

# n = 64
# as = []
# ks = [] 
# for k=2:2:(n+1)
# 	println(k)
# 	# a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	# aftertransition = findall(x->x>0.5, e)
# 	# ac = a[aftertransition[1]]
# 	# println(fmax, 1.0/k)
# 	ac, fmax, iter = cascadebisection(1, 5e-5, n, 1, n-1, 0.0, k, 1.0)
# 	push!(as, ac)
# 	push!(ks, k)
# end 
# plot(ks, as, ".-", color=martaRed, lw=3.0, label="\$n=64\$")

# ks1 = kmns 
# plot(ks, ks1./(ks1.-1.0), "k--", lw=3.0, label="\$\\frac{k}{k-1}\$")
# ylabel("\$\\widetilde{\\alpha}\$", rotation=0, labelpad=20)
# xlabel("\$k\$")
# legend()
tight_layout() 
show()


