using .NetworkResilience
# using PyPlot
# using Statistics
# using Distributed
# using DifferentialEquations
# using LinearAlgebra
# using PyPlot
# using Dagger 
using JuliaDB 

println("Starting...")
# loadtable(glob("powerdata/data/*.csv"), output="powerdatabin", chunks=200)
loadtable("powerdata/Power-Networks-LCL-June2015(withAcornGps)v2_1.csv")
println("Loaded")
