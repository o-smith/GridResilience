using .NetworkResilience
using DifferentialEquations 

#Function to wrap f and ft 
ft(u, p, t) = fswing(u, obj1.adjmat, P, I, D, kappa)
tspan = (0.0, 500.0)
prob = ODEProblem(ft, ψ0, tspan) 

f(du,u,p,t) = fsteadystate(u, A, P, kappa)
df(J,u,p,t) = mainjac(u, E, I, kappa) 

ff = ODEFunction(f;jac=df)
prob = ODEProblem(ff,ones(2),(0.0,10.0))



#Function to wrap df 

#Function to solve

#Initialise network

#Time-step at high coupling 

#Feed initial condition into solver with d and df 

#Iterate kappa down until failure 

#Record critical value  

#Do the above for different realisations 