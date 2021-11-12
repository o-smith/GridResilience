using .NetworkResilience
using PyPlot
using Statistics
using Distributed

#Define colours 
martaRed = "#c24c51" 
martaGreen = "#54a666"
martaBlue = "#4c70b0"
martaPurple = "#7f70b0"
martaGold = "#ccb873"
rc("font", family="serif", size=20)

Ptot = 1.0
k = 6
n = 50
qres = 50 
ensemblesize = 2500 
tol = 5e-4 

ac_corners = zeros(Float64, qres)
ac_middles = zeros(Float64, qres)
qvals = collect(range(0, length=qres, stop=1.0))
fname = "data/qtrans/k6_2500ens.txt"
open(fname, "w") do f 
	for (ind, q) ∈ enumerate(qvals)
		println(q) 
		ac_corner, fmax_corner, iter = cascadebisection(ensemblesize, tol, n, 1, n-1, q, k, Ptot)
		ac_middle, fmax_middle, iter = cascadebisection(ensemblesize, tol, n, div(n,2), div(n,2), q, k, Ptot)

		ac_corners[ind] = ac_corner
		ac_middles[ind] = ac_middle
		
		write(f, "$q $ac_corner $ac_middle \n")

	end
end  

plot(qvals, ac_corners.-ac_middles, "mo", ms=5.0)
show() 