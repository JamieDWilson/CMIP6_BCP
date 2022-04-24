

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

#IPSL-INCA Csoft for historical+ssp1260 and piControl, and their anomaly
# !! the branch of historical run starts 20 years after 1850 of the piControl
download_dir='/gws/pw/j05/cop26_hackathons/bristol/project03/input_nc/multimodel/'
data_dir='/badc/cmip6/data/CMIP6/'
output_dir='/gws/pw/j05/cop26_hackathons/bristol/project03/output_nc/CONTROL/'
Scenario_ssp='ssp126'
Scenario_hist='CMIP'
centre_name='IPSL'
model_name='IPSL-CM5A2-INCA'
exp_pi='piControl'
exp_hist='historical'
exp_ssp='ssp126'
real='r1i1p1f1'
grG='gn'
data_o2 = 'o2'
data_S = 'so'
data_T = 'thetao'
data_V = 'thkcello'

# read the area that is constant and you need for this model
Ar=xr.open_mfdataset('/badc/cmip6/data/CMIP6/CMIP/IPSL/IPSL-CM5A2-INCA/piControl/r1i1p1f1/Ofx/areacello/gn/latest/areacello_Ofx_IPSL-CM5A2-INCA_piControl_r1i1p1f1_gn.nc')

## READ files for piControl run
#oxygen
dir_o2 = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_pi+'/'+real+'/Omon/'+data_o2+'/'+grG+'/latest/'
file_o2=glob.glob(dir_o2+'/*.nc')
o2_list=(file_o2[0])
data_pi_o2 = xr.open_mfdataset(o2_list)
#temperature
dir_T = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_pi+'/'+real+'/Omon/'+data_T+'/'+grG+'/latest/'
file_T=glob.glob(dir_T+'/*.nc')
T_list=(file_T[0])
data_pi_T = xr.open_mfdataset(T_list)
#Salinity
data_pi_S = xr.open_mfdataset(download_dir+'piControl/'+data_S+'_Omon_'+model_name+'_'+exp_pi+'_'+real+'_'+grG+'_'+'*.nc')
#Volcello - Volume
data_pi_V = xr.open_mfdataset(download_dir+'piControl/'+data_V+'_Omon_'+model_name+'_'+exp_pi+'_'+real+'_'+grG+'_'+'*.nc')
data_pi_V=fs.area_to_volume(Ar,data_pi_V)
data_pi_V=data_pi_V.transpose("time", "deptht", "y","x")

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
#Volcello - Volume
dir_V = data_dir+Scenario_hist+'/'+centre_name+'/'+model_name+'/'+exp_hist+'/'+real+'/Omon/'+data_V+'/'+grG+'/latest/'
data_hist_V = xr.open_mfdataset(dir_V+'/*.nc')
data_hist_V=fs.area_to_volume(Ar,data_hist_V)
data_hist_V=data_hist_V.transpose("time", "deptht", "y","x")

## READ files for ssp126 run
#oxygen
data_ssp_o2 = xr.open_mfdataset(download_dir+Scenario_ssp+'/'+data_o2+'_Omon_'+model_name+'_'+exp_ssp+'_'+real+'_'+grG+'_'+'*.nc')
#temperature
data_ssp_T = xr.open_mfdataset(download_dir+Scenario_ssp+'/'+data_T+'_Omon_'+model_name+'_'+exp_ssp+'_'+real+'_'+grG+'_'+'*.nc')
#Salinity
data_ssp_S = xr.open_mfdataset(download_dir+Scenario_ssp+'/'+data_S+'_Omon_'+model_name+'_'+exp_ssp+'_'+real+'_'+grG+'_'+'*.nc')
#Volcello - Volume
file_V=glob.glob(download_dir+Scenario_ssp+'/'+data_V+'_Omon_'+model_name+'_'+exp_ssp+'_'+real+'_'+grG+'_'+'*.nc')
data_ssp_V = xr.open_mfdataset(file_V)
data_ssp_V=fs.area_to_volume(Ar,data_ssp_V)
data_ssp_V=data_ssp_V.transpose("time", "deptht", "y","x")

## Find the time for historical+ssp run and set accordingly nstart and end for piControl
# !! Problem with IPSL-INCA pi control not all I need (only 330 years of actuall picontrol that I can use 1870-2099)
# so INCA will be estimated only up to 2080 rather than 2100 
time_all=xr.concat([data_hist_o2.time,data_ssp_o2.time],dim='time')
nstart1=240 # IPSL-INCA starts at 20 years after 1850, so 240 months (timesteps) after the first file
nend1=time_all.size+nstart1-21*12

## Caclulate Csoft control historical and ssp
Csoft_pi=fs.compute_AOU_Csoft(data_pi_o2,data_pi_T,data_pi_S,data_pi_V,model_name,nstart1,nend1,output_dir)

nstart=0
nend2=data_hist_o2.time.size
Csoft_hist=fs.compute_AOU_Csoft(data_hist_o2,data_hist_T,data_hist_S,data_hist_V,model_name,nstart,nend2,output_dir)

nstart=0
nend3=data_ssp_o2.time.size-21*12
Csoft_ssp=fs.compute_AOU_Csoft(data_ssp_o2,data_ssp_T,data_ssp_S,data_ssp_V,model_name,nstart,nend3,output_dir)

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

Csoft_all.to_netcdf(output_dir+'Csoft_yearly_all_'+model_name+'_ssp126.nc',mode='w')
