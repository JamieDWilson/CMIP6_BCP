

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

#CMCC-ESM2 Csoft for historical+ssp370 and piControl, and their anomaly
# !! the branch of historical run starts from year 1850 of the piControl
download_dir='/gws/pw/j05/cop26_hackathons/bristol/project03/input_nc/multimodel/'
data_dir='/badc/cmip6/data/CMIP6/'
output_dir='/gws/pw/j05/cop26_hackathons/bristol/project03/output_nc/CONTROL/'
Scenario_ssp='ScenarioMIP'
Scenario_hist='CMIP'
centre_name='CMCC'
model_name='CMCC-ESM2'
exp_pi='piControl'
exp_hist='historical'
exp_ssp='ssp370'
real='r1i1p1f1'
grG='gn'
data_o2 = 'o2'
data_S = 'so'
data_T = 'thetao'
data_V = 'volcello'


## READ files for piControl run
#oxygen
dir_o2 = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_pi+'/'+real+'/Omon/'+data_o2+'/'+grG+'/latest/'
file_o2=glob.glob(dir_o2+'/*.nc')
o2_list=(file_o2[0:13])
data_pi_o2 = xr.open_mfdataset(o2_list)
#temperature
dir_T = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_pi+'/'+real+'/Omon/'+data_T+'/'+grG+'/latest/'
file_T=glob.glob(dir_T+'/*.nc')
T_list=(file_T[0:13])
data_pi_T = xr.open_mfdataset(T_list)
#Salinity
data_pi_S = xr.open_mfdataset(download_dir+'piControl/'+data_S+'_Omon_'+model_name+'_'+exp_pi+'_'+real+'_'+grG+'_'+'*.nc')

## READ files for historical run
#oxygen
dir_o2 = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_hist+'/'+real+'/Omon/'+data_o2+'/'+grG+'/latest/'
data_hist_o2 = xr.open_mfdataset(dir_o2+'/*.nc')
#temperature
dir_T = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_hist+'/'+real+'/Omon/'+data_T+'/'+grG+'/latest/'
data_hist_T = xr.open_mfdataset(dir_T+'/*.nc')
#Salinity
dir_S = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_hist+'/'+real+'/Omon/'+data_S+'/'+grG+'/latest/'
data_hist_S = xr.open_mfdataset(dir_S+'/*.nc')


## READ files for ssp370 run
#oxygen
dir_o2 = data_dir+Scenario_ssp+'/'+centre_name+'/'+model_name+'/'+exp_ssp+'/'+real+'/Omon/'+data_o2+'/'+grG+'/latest/'
data_ssp_o2 = xr.open_mfdataset(dir_o2+'/*.nc')
#temperature
dir_T = data_dir+Scenario_ssp+'/'+centre_name+'/'+model_name+'/'+exp_ssp+'/'+real+'/Omon/'+data_T+'/'+grG+'/latest/'
data_ssp_T = xr.open_mfdataset(dir_T+'/*.nc')
#Salinity
dir_S = data_dir+Scenario_ssp+'/'+centre_name+'/'+model_name+'/'+exp_ssp+'/'+real+'/Omon/'+data_S+'/'+grG+'/latest/'
data_ssp_S = xr.open_mfdataset(dir_S+'/*.nc')

## Volcello - Volume, single file as a constalnt volume
data_V = xr.open_mfdataset('/badc/cmip6/data/CMIP6/CMIP/CMCC/CMCC-ESM2/piControl/r1i1p1f1/Ofx/volcello/gn/latest/volcello_Ofx_CMCC-ESM2_piControl_r1i1p1f1_gn.nc')

## Find the time for historical+ssp run and set accordingly nstart and end for piControl
time_all=xr.concat([data_hist_o2.time,data_ssp_o2.time],dim='time')
nstart=0 # CMCC-ESM2  starts at year 1850,at the beginning of file 1850 which is the first file read
nend1=time_all.size

## Caclulate Csoft control historical and ssp
Csoft_pi=fs.compute_AOU_Csoft(data_pi_o2,data_pi_T,data_pi_S,data_V,model_name,nstart,nend1,output_dir)

nend2=data_hist_o2.time.size
Csoft_hist=fs.compute_AOU_Csoft(data_hist_o2,data_hist_T,data_hist_S,data_V,model_name,nstart,nend2,output_dir)

nend3=data_ssp_o2.time.size
Csoft_ssp=fs.compute_AOU_Csoft(data_ssp_o2,data_ssp_T,data_ssp_S,data_V,model_name,nstart,nend3,output_dir)

## merge Csoft historical and ssp into one timeserie
Csoft_cl=xr.concat([Csoft_hist,Csoft_ssp],dim='year')

#Csoft_control=xr.DataArray(Csoft_pi,name='Csoft_control',coords=[time_all[:]], dims=["time"])
## Write output in netcdf file
Csoft_all = xr.Dataset(
    {
        "Csoft_control": (["year"], Csoft_pi.Csoft),
        "Csoft_climate": (["year"], Csoft_cl.Csoft),
    },
    coords={
        "year" : Csoft_cl.year,
    },
)

Csoft_all.to_netcdf(output_dir+'Csoft_yearly_all_'+model_name+'.nc',mode='w')
