

#!/usr/bin/env python3
#Control run before this script: module load jaspy 
# run before this script: source /gws/pw/j05/cop26_hackathons/bristol/activate-env

from itertools import chain
import glob

import os
import xarray as xr
import numpy as np
import seawater as sw
from netCDF4 import Dataset,num2date,date2num

import hackathon_functions_yearly as fs

#GFDL-ESM4 Csoft for historical+ssp370 and piControl, and their anomaly
# !! the branch of historical run starts from year 101 of the piControl
data_dir='/badc/cmip6/data/CMIP6/'
output_dir='/gws/pw/j05/cop26_hackathons/bristol/project03/output_nc/CONTROL/'
Scenario_ssp='ScenarioMIP'
Scenario_hist='CMIP'
centre_name='NOAA-GFDL'
model_name='GFDL-ESM4'
exp_pi='piControl'
exp_hist='historical'
exp_ssp='ssp370'
real='r1i1p1f1'
grG='gr'
data_o2 = 'o2'
data_S = 'so'
data_T = 'thetao'

Ar=xr.open_mfdataset('/badc/cmip6/data/CMIP6/CMIP/NOAA-GFDL/GFDL-ESM4/piControl/r1i1p1f1/Ofx/areacello/gr/latest/areacello_Ofx_GFDL-ESM4_piControl_r1i1p1f1_gr.nc')

## READ files for piControl run
#temperature
dir_T = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_pi+'/'+real+'/Omon/'+data_T+'/'+grG+'/latest/'
file_T=glob.glob(dir_T+'/*.nc')
T_list=(file_T[5:19])
data_pi_T = xr.open_mfdataset(T_list)
#Salinity
dir_S = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_pi+'/'+real+'/Omon/'+data_S+'/'+grG+'/latest/'
file_S=glob.glob(dir_S+'/*.nc')
S_list=(file_S[5:19])
data_pi_S = xr.open_mfdataset(S_list)

## READ files for historical run
#temperature
dir_T = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_hist+'/'+real+'/Omon/'+data_T+'/'+grG+'/latest/'
data_hist_T = xr.open_mfdataset(dir_T+'/*.nc')
#Salinity
dir_S = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_hist+'/'+real+'/Omon/'+data_S+'/'+grG+'/latest/'
data_hist_S = xr.open_mfdataset(dir_S+'/*.nc')

## READ files for ssp370 run
#temperature
dir_T = data_dir+Scenario_ssp+'/'+centre_name+'/'+model_name+'/'+exp_ssp+'/'+real+'/Omon/'+data_T+'/'+grG+'/latest/'
data_ssp_T = xr.open_mfdataset(dir_T+'/*.nc')
#Salinity
dir_S = data_dir+Scenario_ssp+'/'+centre_name+'/'+model_name+'/'+exp_ssp+'/'+real+'/Omon/'+data_S+'/'+grG+'/latest/'
data_ssp_S = xr.open_mfdataset(dir_S+'/*.nc')


## Find the time for historical+ssp run and set accordingly nstart and end for piControl
time_all=xr.concat([data_hist_T.time,data_ssp_T.time],dim='time')
nstart=0 # GFDL starts at year 101,at the beginning of file 101 which is the first file read
nend1=time_all.size

## Caclulate Csoft control historical and ssp
Csoft_pi=fs.compute_stratification_timeseries(data_pi_T,data_pi_S,Ar,model_name,nstart,nend1,output_dir)

nend2=data_hist_T.time.size
Csoft_hist=fs.compute_stratification_timeseries(data_hist_T,data_hist_S,Ar,model_name,nstart,nend2,output_dir)

nend3=data_ssp_T.time.size
Csoft_ssp=fs.compute_stratification_timeseries(data_ssp_T,data_ssp_S,Ar,model_name,nstart,nend3,output_dir)

## merge Csoft historical and ssp into one timeserie
#Csoft_cl=xr.concat([Csoft_hist.Csoft,Csoft_ssp.Csoft],dim='year')
Csoft_cl=xr.concat([Csoft_hist,Csoft_ssp],dim='year')
#Csoft_cl=Csoft_cl.to_dataset()

#Csoft_control=xr.DataArray(Csoft_pi,name='Csoft_control',coords=[time_all[:]], dims=["time"])
## Write output in netcdf file
#time_years=np.arange(1850,1850+nend1/12,1/12)
#time_new=time_all.groupby('time.year').mean('time').compute()
Csoft_all = xr.Dataset(
    {
        "Csoft_control": (["year"], Csoft_pi.stratification),
        "Csoft_climate": (["year"], Csoft_cl.stratification),
    },
    coords={
        "year" : Csoft_cl.year,
    },
)

Csoft_all.to_netcdf(output_dir+'Strat_yearly_all_'+model_name+'_ssp370.nc',mode='w')
