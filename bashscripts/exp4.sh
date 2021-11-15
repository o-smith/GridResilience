#!/bin/bash
echo Running experiment to compute cascade durations at different parameter values...
mkdir -p data/durations
export OMP_NUM_THREADS=1
julia -L src/NetworkResilience.jl scripts/responsetimes.jl 
python scripts/responsesplot.py 