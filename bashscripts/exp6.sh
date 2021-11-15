#!/bin/bash
echo Running script to \read \in the whole power consumption and PV data sets and plot their mean timeseries... 
#Check that the power data has been downloaded
if [ -d "powerdata/PV" ] && [ -d "powerdata/data" ]; then
  python scripts/timeseriesplot.py
else
  # No power data has been downloaded
  echo Warning: power data has not been downloaded.
  echo Please download PV data from:
  echo https://data.london.gov.uk/dataset/photovoltaic--pv--solar-panel-energy-generation-data  
  echo and place it \in a directory called \"powerdata/PV\".  
  echo Please download the power consumption data from: 
  echo https://data.london.gov.uk/dataset/smartmeter-energy-use-data-in-london-households 
  echo making sure to \select \the 168 file option, and \then place it \in a directory called \"powerdata/data\".  
fi
