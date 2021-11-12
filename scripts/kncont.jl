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
f = [] 
ns = []
k = 4
for n=6:2:200
	println(n)
	ac, fmax, iter = cascadebisection(1, 5e-5, n, 1, n-1, 0.0, k, 1.0)
	push!(f, fmax)
	push!(ns, n)
end 
plot(f)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
# plot((1.0/k).*x)
x = (ns.-k)./(ns.^2)
plot(0.25.*(1.0.+x))
tight_layout() 
show() 
# n = 50
# as = []
# ks = [] 
# for k=4:2:(n-4)
# 	println(k)
# 	# a, e, std, t, tstd, mm, mato, fmax = smcascadefixedtop(1, n, 1, n-1, 0.0, 4000, 0.01, 10.0, k, 1.0, 1.0)
# 	# aftertransition = findall(x->x>0.5, e)
# 	# ac = a[aftertransition[1]]
# 	# println(fmax, 1.0/k)
# 	ac, fmax, iter = cascadebisection(1, 5e-5, n, 1, n-1, 0.0, k, 1.0)
# 	push!(as, fmax)
# 	push!(ks, k)
# end
# ks1 = ks  
# plot(ks, as.-(1.0./ks1), ".-", color=martaGreen, lw=3.0, label="\$n=16\$")


# # plot(ks1, 1.0./(10.0.*ks1), "k--", lw=3.0, label="\$\\frac{1}{k}\$")
# ylabel("\$\\widetilde{\\alpha}\$", rotation=0, labelpad=20)
# xlabel("\$k\$")
# legend()
# tight_layout()
# show() 