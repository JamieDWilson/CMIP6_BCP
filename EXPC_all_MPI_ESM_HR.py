

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

#MPI-ESM-HR expc for historical+ssp370 and piControl, and their anomaly
# !! the branch of historical run starts from year 1850 of the piControl
download_dir='/gws/pw/j05/cop26_hackathons/bristol/project03/input_nc/multimodel/'
data_dir='/badc/cmip6/data/CMIP6/'
output_dir='/gws/pw/j05/cop26_hackathons/bristol/project03/output_nc/CONTROL/'
Scenario_ssp='ScenarioMIP'
Scenario_hist='CMIP'
centre_name='MPI-M'
model_name='MPI-ESM1-2-HR'
exp_pi='piControl'
exp_hist='historical'
exp_ssp='ssp370'
real='r1i1p1f1'
grG='gn'
data_expc = 'expc'

# Give here the year you want the map. In this case the pre-industrial is actually from the piControl not 1850-1860
#timeslices=[2050,2060]
timeslices=[2090,2100]

Ar=xr.open_mfdataset('/badc/cmip6/data/CMIP6/CMIP/MPI-M/MPI-ESM1-2-HR/piControl/r1i1p1f1/Ofx/areacello/gn/latest/areacello_Ofx_MPI-ESM1-2-HR_piControl_r1i1p1f1_gn.nc')
## READ files for piControl run
#expc
file_pi_expc=glob.glob(download_dir+'piControl/'+data_expc+'_Omon_'+model_name+'_'+exp_pi+'_'+real+'_'+grG+'_'+'*.nc')
expc_list=(file_pi_expc[0:51])
data_pi_expc = xr.open_mfdataset(expc_list)

## READ files for historical run
#expc
file_hist_expc=glob.glob(download_dir+'historical/'+data_expc+'_Omon_'+model_name+'_'+exp_hist+'_'+real+'_'+grG+'_'+'*.nc')
data_hist_expc = xr.open_mfdataset(file_hist_expc)

## READ files for ssp370 run
#expc
file_ssp_expc=glob.glob(download_dir+'ssp370/'+data_expc+'_Omon_'+model_name+'_'+exp_ssp+'_'+real+'_'+grG+'_'+'*.nc')
data_ssp_expc = xr.open_mfdataset(file_ssp_expc)

## Find the time for historical+ssp run and set accordingly nstart and end for piControl
time_all=xr.concat([data_hist_expc.time,data_ssp_expc.time],dim='time')
nstart=0 # MPI starts at year 1850,at the beginning of file 1850 which is the first file read
nend1=time_all.size
# MPI_LR 1850 is year 1850 
datch=0
timeslices_pi=[element-datch for element in timeslices]

## Caclulate expc control historical and ssp
EXPC_pi=fs.compute_transfer_efficiency(data_pi_expc,Ar,timeslices_pi,datch,model_name,exp_pi,nstart,nend1,output_dir)

nend2=data_hist_expc.time.size
timeslices_hist=[1850,1860]#just go with a historical start
EXPC_hist=fs.compute_transfer_efficiency(data_hist_expc,Ar,timeslices_hist,0,model_name,exp_hist,nstart,nend2,output_dir)

nend3=data_ssp_expc.time.size
timeslices_ssp=timeslices#if it is within the ssp period
EXPC_ssp=fs.compute_transfer_efficiency(data_ssp_expc,Ar,timeslices_ssp,0,model_name,exp_ssp,nstart,nend3,output_dir)

## merge expc historical and ssp into one timeserie
EXPC_cl=xr.concat([EXPC_hist[0],EXPC_ssp[0]],dim='year')
EXPC_control=EXPC_pi[0]

## Write output in netcdf file
POC_all = xr.Dataset(
    {
        "expc_control": (["year"], EXPC_control.expc),
        "expc_climate": (["year"], EXPC_cl.expc),
        "exp_control": (["year"], EXPC_control.exp),
        "exp_climate": (["year"], EXPC_cl.exp),
        "deep_control": (["year"], EXPC_control.deep),
        "deep_climate": (["year"], EXPC_cl.deep),
    },
    coords={
        "year" : EXPC_cl.year,
    },
)

POC_all.to_netcdf(output_dir+'POC_yearly_all_'+model_name+'.nc',mode='w')

