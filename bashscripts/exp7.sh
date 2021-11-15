#!/bin/bash
export OMP_NUM_THREADS=1 
echo Running script to generate ensembles of micro-grids from data sampled from the PV and consumption data sets,
echo this will construct weekly demand-generation profiles \for different levels of uptake, 
echo different \times of year, and with and without batteries. 
#Check that the power data has been downloaded
if [ -d "powerdata/PV" ] && [ -d "powerdata/data" ]; then
  echo First, generating the ensembles from the data 
  python scripts/powerexperiments.py
else
  echo Warning: power data has not been downloaded.
  echo Please download PV data from:
  echo https://data.london.gov.uk/dataset/photovoltaic--pv--solar-panel-energy-generation-data  
  echo and place it \in a directory called \"powerdata/PV\".  
  echo Please download the power consumption data from: 
  echo https://data.london.gov.uk/dataset/smartmeter-energy-use-data-in-london-households 
  echo making sure to \select \the 168 file option, and \then place it \in a directory called \"powerdata/data\".  
  python scripts/powerexperiments.py
fi 