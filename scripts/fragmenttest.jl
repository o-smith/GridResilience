using .NetworkResilience
using Plots
using Statistics
using Distributed

alpha = 0.8
ensenmblesize = 500
fragvec = []
anim = @animate for alpha=0.7:0.05:1.1
	println(alpha)
	for i=1:ensenmblesize
		#println(i)
		frags = simplecascade(100, 50, 50, 0.1, alpha, 4, 1.0, 1.0)
		push!(fragvec, frags)
	end
	histogram(fragvec, bins=50)
	#xrange(0,200)
	#yrange(0,250)
	#xlabel("Fragments")
end

gif(anim, "q1frags_50_50_0.gif", fps = 15)
# show() 
