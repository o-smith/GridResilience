using .NetworkResilience

ensemblesize = 200
tol = 5e-4 
n = 60
ns = 15
nd = 45
q = 0.1
I = 1.0 
D = 1.0 
k = 4 
Ptot = 1.0 
κ = 5.0 

#Do cascade as a function of capacity alpha and 
#save the fractions of nodes failing by overload or desync
println("Computing cascades for case (n+,n-,np)=(15,45,0)")
alphas, s, ov, desy, deb, fmaxes, crit = swingcascadebreakdownvsalpha(ensemblesize, n, ns, nd, q, k, 50, 0.1, 2.5, I, D, κ)  
fname = "data/cascadefailurebreakdown1.txt"
open(fname, "w") do f 
    for i=1:length(alphas)
        beta = alphas[i]
        gamma = s[i]
        kappa = ov[i]
        epsil = desy[i] + deb[i]
        write(f, "$beta $gamma $kappa $epsil \n")
    end
end

#And again, but this time for the perfectly distributed generation case
println("Computing cascades for case (n+,n-,np)=(30,30,0)")
ns = 30
nd = 30
alphas, s, ov, desy, deb, fmaxes, crit = swingcascadebreakdownvsalpha(ensemblesize, n, ns, nd, q, k, 50, 0.1, 2.5, I, D, κ)  
fname = "data/cascadefailurebreakdown2.txt"
open(fname, "w") do f 
    for i=1:length(alphas)
        beta = alphas[i]
        gamma = s[i]
        kappa = ov[i]
        epsil = desy[i] + deb[i]
        write(f, "$beta $gamma $kappa $epsil \n")
    end
end

#Now compute mean cascade duration as a function of alpha
println("Computing cascade durations for case (n+,n-,np)=(15,45,0)")
ns = 15
nd = 45
a, s, t, fmax, crits = swingcascadeandtvalpha(ensemblesize, n, ns, nd, q, k, 50, 0.1, 2.5, I, D, κ)
fname = "data/cascadedurations1.txt"
open(fname, "w") do f 
    for i=1:length(a)
        beta = a[i]
        gamma = t[i]
        write(f, "$beta $gamma \n")
    end
end

#And again, but this time for the perfectly distributed generation case
println("Computing cascade durations for case (n+,n-,np)=(30,30,0)")
ns = 15
nd = 45
a, s, t, fmax, crits = swingcascadeandtvalpha(ensemblesize, n, ns, nd, q, k, 50, 0.1, 2.5, I, D, κ)
fname = "data/cascadedurations2.txt"
open(fname, "w") do f 
    for i=1:length(a)
        beta = a[i]
        gamma = t[i]
        write(f, "$beta $gamma \n")
    end
end


