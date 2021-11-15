using .NetworkResilience

ensemblesize = 500 
n = 50 
q = 0.1 
ns = 10 
nd = 40 
k = 4
I = 1.0 
D = 1.0 
kappa = 5.0  

println("Computing cascades for alpha=0.75")
alpha = 0.75
s, t = svt(ensemblesize, n, ns, nd, q, k, alpha, I, D, kappa)
fname = "data/cascade_s_vs_t_q75.txt"
open(fname, "w") do f 
    for i=1:length(s)
        beta = s[i]
        gamma = t[i]
        write(f, "$beta $gamma \n")
    end
end

println("Computing cascades for alpha=1")
alpha = 1.0
s, t = svt(ensemblesize, n, ns, nd, q, k, alpha, I, D, kappa)
fname = "data/cascade_s_vs_t_q10.txt"
open(fname, "w") do f 
    for i=1:length(s)
        beta = s[i]
        gamma = t[i]
        write(f, "$beta $gamma \n")
    end
end

println("Computing cascades for alpha=1.21")
alpha = 1.21
s, t = svt(ensemblesize, n, ns, nd, q, k, alpha, I, D, kappa)
fname = "data/cascade_s_vs_t_q121.txt"
open(fname, "w") do f 
    for i=1:length(s)
        beta = s[i]
        gamma = t[i]
        write(f, "$beta $gamma \n")
    end
end

println("Computing cascades for alpha=2")
alpha = 2.0
s, t = svt(ensemblesize, n, ns, nd, q, k, alpha, I, D, kappa)
fname = "data/cascade_s_vs_t_q20.txt"
open(fname, "w") do f 
    for i=1:length(s)
        beta = s[i]
        gamma = t[i]
        write(f, "$beta $gamma \n")
    end
end