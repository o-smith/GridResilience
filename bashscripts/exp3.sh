#!/bin/bash 
echo Running expermiment to compute distributions of rho...
export OMP_NUM_THREADS=1
julia -L src/NetworkResilience.jl scripts/swingdists.jl 
python scripts/distsplot.py