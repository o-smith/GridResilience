# Code in support of the paper "The effect of renewable energy incorporation on power grid stability and resilience" 


[![DOI](https://zenodo.org/badge/427405018.svg)](https://zenodo.org/badge/latestdoi/427405018)



Respository containing data and code for "The effect of renewable energy incorporation on power grid stability and resilience" 

## Requirements

The code requires Python with the packages: 

- numpy
- scipy 
- matplotlib 
- python-ternary 
- dill
- dask 
- pandas 

All of the above can be installed using pip, e.g.: 

```bash
pip install dill 
```

Please also run the update script:
```bash
python -m pip install "dask[dataframe]" --upgrade 
```

The code also requires Julia, version 1 or later with packages: 

- LinearAlgebra 
- Statistics 
- Distributed 
- Distributions 
- DifferentialEquations 

---
## Usage 

Bash scripts have been included in the bashscripts directory to run simulations and produce plots from the paper. To run them, first ensure they are executable; e.g.
```bash
chmod +x bashscripts/exp1.sh
```

The scripts default to using pre-computed micro-grid trajectory data. To run the simulations to generate the trajctory data from scratch from the household consumption data sets, please download the household consumption and PV generation data from:

https://data.london.gov.uk/dataset/smartmeter-energy-use-data-in-london-households 

and 

https://data.london.gov.uk/dataset/photovoltaic--pv--solar-panel-energy-generation-data 


respectively. These should be placed in directories called powerdata/data powerdata/PV directory respectively. The scripts also default to using pre-computed simplex data. To run the simulations to generate the simplex data from scratch, run the bash scripts stored in the hpcscsripts directory warning: the bash scripts within the hpcscripts folder are computationally demanding and are inteded for use on an HPC.)

---

Contact oliver.smith@nottingham.ac.uk for any further assistance running scripts. 
