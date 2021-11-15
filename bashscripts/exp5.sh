#!/bin/bash
echo Running experiment to compute cascade durations at different alpha values...
export OMP_NUM_THREADS=1
julia -L src/NetworkResilience.jl scripts/timeexperiments.jl
python scripts/scatterhists.py
