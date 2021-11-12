using .NetworkResilience
using PyPlot
using Statistics
using Distributed

n = 200
z = 100

alphamin = 0.01
alphamax = 3.0
alphares = 300
as = zeros(Float64, z-1, alphares)
es = zeros(Float64, z-1, alphares)
# for i=2:z
# 	a, e, std, t, tstd, mat = smcascadeclock(n, i, 0.0, alphares, alphamin, alphamax, 2, 1.0, 1.0)
# 	fill_between(a, e.+std, e.-std, color="k", alpha=0.2)
# 	as[i-1,:] = a 
# 	es[end-i+2,:] = e
# 	println(i)
# end 

for i=2:z
	e = smcascadeclockrand(n, i, 0.0, alphares, alphamin, alphamax, 2, 1.0, 1.0)
	es[end-i+2,:] = e
	println(i)
end 
# fname = "data/rho_dostancek10_200_normed.txt"
# open(fname, "w") do f 
# 	for i=1:z-1
# 		for j=1:alphares
# 			beta = es[i,j]
# 			write(f, "$beta ")
# 		end 
# 		write(f, "\n")
# 	end
# end

x = collect(range(0.01, length=alphares, stop=1.0)) 
# s = 62.0.*x .+ 11 .- 31
# plot(x, 64.0.*es[20,:])
# plot(x, s)
# show() 


stepsize = (alphamax - alphamin)/alphares
imshow(es, interpolation="none", cmap="jet", aspect="auto", extent=[alphamin,alphamax,0,z-2]) #, vmin=0, vmax=1)
# plot([2,2], [10,40], "k")

# curve = floor.(x.*(n-2)) .+ 0.2  
# curve2 = floor.(x.*(2-n) .+ (n)/2) .- 0.5  
# plot(x, curve,  "w-", lw=2.0)
# plot(x, curve2,  "w-", lw=2.0)
xlim([alphamin,alphamax])
ylim([0,z-2])
xlabel("Alpha", fontsize=18)
ylabel("Source separation distance", fontsize=18)
colorbar() 
show() 

