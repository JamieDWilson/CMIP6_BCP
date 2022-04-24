# CMIP6_BCP

This repository containts the Python code to generate the analysis for Wilson et al., "The Biological Carbon Pump in CMIP6 models: 21st century trends and uncertainties".

1) hackathon_functions_yearly.py:
  - individual analysis functions (see also comments within file)
    - dim_info: converts model specific dimension names to common names
    - compute_AOU_Csoft: calculates annual mean volume-integrated Apparent Oxygen Utilisation timeseries 
    - compute_o2sat: calculates oxygen saturation
    - check_inputs: checks and converts dimension units
    - area_to_volume: converts area to volume when explicit output not available
    - compute_transfer_efficiency: calculates annual mean timeseries of POC fluxes at 100m, 1000m and transfer efficiency, plus time-averaged spatial fields
  
2) EXPC_all_*.py:
  - specific scripts to calculate outputs related to POC fluxes from an individual model
  
3) Csoft_all_*.py: 
  - specific scripts to calculate outputs related to Apparent Oxygen Utilisation and Csoft from an individual model

4) Strat_all_*.py
  - specific scripts to calculate outputs related to stratification from an individual model
  
 
Authors: Jamie D. Wilson (University of Bristol) & Anna Katavouta (National Oceanography Centre) adapting code from all authors.
  
