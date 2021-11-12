using .NetworkResilience
using PyPlot
using Statistics

# n = 64
# k = 4
# q = 1.0
# tol = 5e-5
# ensemblesize = 50 
# m = div(n*k,2)
# rvec = zeros(Float64, m*ensemblesize) 
# avec = zeros(Float64, m*ensemblesize) 
# for i=1:ensemblesize
# 	println(i)
# 	ns = 32 #rand(1:div(n,2))
# 	nd = n - ns
# 	r, a, fmax = sequencedstartbisection(tol, n, ns, nd, q, k) 
# 	rvec[((i-1)*m)+1:i*m] = r 
# 	avec[((i-1)*m)+1:i*m] = a 
# end 


# plot(rvec, avec, "ko") 
# show() 
# hist(avec, bins=50)
# show() 
# x = cor(r, a)
# println(x)

# hist(a)
# show() 
avec = zeros(Float64, 2000)

fname = "data/achists/sf_30_30_m3_cascade2.txt"
open(fname, "w") do f 

	for i=1:2000
		println(i)
		a, fmax, iter = cascadebisection2(1, 5e-5, 60, 30, 30, 0.1, 4, 1.0)
		avec[i] = a 
		write(f, "$a \n")
	end
	
end

println(mean(avec))
hist(avec, bins=20)
show() 
