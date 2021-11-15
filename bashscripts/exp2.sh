#!/bin/bash
echo Running experiment to analyse cascade failure breakdown modes...
export OMP_NUM_THREADS=1
julia -L src/NetworkResilience.jl scripts/swingdriver.jl 
python scripts/tempplot.py 