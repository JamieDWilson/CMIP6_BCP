

#!/usr/bin/env python3
# run before this script: module load jaspy 
# run before this script: source /gws/pw/j05/cop26_hackathons/bristol/activate-env

from itertools import chain
import glob

import os
import xarray as xr
import numpy as np

import hackathon_functions as hfs

data_dir='../../input_nc/multimodel/historical/'
output_dir='../../output_nc/master/o2/'
mip='ScnenarioMIP'
data='o2'
exp='historical'
timeslices=['1850','1860']
#timeslices=['2015','2025','2090','2100']


# CESM2-FVR
# -> OK! 
#centre_name='NCAR'
#model_name='CESM2-FV2'
#data_list=glob.glob(data_dir+'expc_*_'+model_name+'_*')
#grid_list=['/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/'+exp+'/r1i1p1f1/Ofx/areacello/gn/latest/areacello_Ofx_'+model_name+'_'+exp'_r1i1p1f1_gn.nc']
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir)

# CESM2
# ->  OK!
#centre_name='NCAR'
#model_name='CESM2'
#data_list=glob.glob(data_dir+'expc_*_'+model_name+'_*')
#grid_list=['/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Ofx/areacello/gr/latest/areacello_Ofx_'+model_name+'_historical_r1i1p1f1_gr.nc']
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir)

# CESM2-WACCM-FV2
# -> OK!
#centre_name='NCAR'
#model_name='CESM2-WACCM-FV2'
#data_list=glob.glob(data_dir+'expc_*_'+model_name+'_*')
#grid_list=['/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Ofx/areacello/gn/latest/areacello_Ofx_'+model_name+'_historical_r1i1p1f1_gn.nc']
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir)

# CESM2-WACCM
# -> missing areacello
#centre_name='NCAR'
#model_name='CESM2-WACCM'
#data_list=glob.glob(data_dir+'expc_*_'+model_name+'_*')
#grid_list=[data_dir+'areacello_Ofx_'+model_name+'_historical_r1i1p1f1_gn.nc']
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir)

# GFDL-CM4
# -> missing areacello 
#centre_name='NOAA-GFDL'
#model_name='GFDL-CM4'
#data_list=glob.glob(data_dir+'o2_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so_*_'+model_name+'_*')
#grid_list=glob.glob(data_dir+'volcello_*_'+model_name+'_*')
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list)

# GFDL-ESM4
# -> OK!
# -> OK! (using historical areacello file - no ssp370 in badc)
#centre_name='NOAA-GFDL'
#model_name='GFDL-ESM4'
#data_list=glob.glob(data_dir+data+'_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao'+'_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so'+'_*_'+model_name+'_*')
#grid_list=glob.glob(data_dir+'volcello'+'_*_'+model_name+'_*')
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list)

# IPSL-CM5A2-INCA
# -> OK!
# -> OK! 
#centre_name='IPSL'
#model_name='IPSL-CM5A2-INCA'
#data_list=glob.glob(data_dir+data+'_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao'+'_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so'+'_*_'+model_name+'_*')
#grid_list=['/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Ofx/areacello/gn/latest/areacello_Ofx_'+model_name+'_historical_r1i1p1f1_gn.nc']
#thk_list=glob.glob('/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Omon/thkcello/gn/latest/*')
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list,thk_list)

# IPSL-CM6A-LR
# -> 
# -> OK!
#centre_name='IPSL'
#model_name='IPSL-CM6A-LR'
#data_list=glob.glob(data_dir+data+'_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao'+'_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so'+'_*_'+model_name+'_*')
#grid_list=['/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Ofx/areacello/gn/latest/areacello_Ofx_'+model_name+'_historical_r1i1p1f1_gn.nc']
#thk_list=glob.glob('/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Omon/thkcello/gn/latest/*')
#print(thk_list)
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list,thk_list)

# IPSL-CM6A-LR-INCA
# -> just seems to hang when interpolating?! 
centre_name='IPSL'
model_name='IPSL-CM6A-LR-INCA'
data_list=glob.glob(data_dir+'o2_*_'+model_name+'_*')
t_list=glob.glob(data_dir+'thetao_*_'+model_name+'_*')
s_list=glob.glob(data_dir+'so_*_'+model_name+'_*')
grid_list=glob.glob(data_dir+'areacello_*_'+model_name+'_*')
thk_list=glob.glob(data_dir+'thkcello_*_'+model_name+'_*')
hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list,thk_list)

# MPI-ESM-1-2-HAM
# -> OK!
# -> OK!
centre_name='HAMMOZ-Consortium'
#model_name='MPI-ESM-1-2-HAM'
#data_list=glob.glob(data_dir+'o2_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so_*_'+model_name+'_*')
#grid_list=glob.glob(data_dir+'volcello'+'_*_'+model_name+'_*')
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list)

# MPI-EXM1_2_HR
# -> OK!
# -> OK! 
#centre_name='MPI-M'
#model_name='MPI-ESM1-2-HR'
#data_list=glob.glob(data_dir+'o2_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so_*_'+model_name+'_*')
#grid_list=['/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Ofx/volcello/gn/latest/volcello_Ofx_'+model_name+'_historical_r1i1p1f1_gn.nc']
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list)

# MPI-EXM1_2_LR
# -> OK! 
#centre_name='MPI-M'
#model_name='MPI-ESM1-2-LR'
#data_list=glob.glob(data_dir+'o2_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so_*_'+model_name+'_*')
#grid_list=['/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Ofx/volcello/gn/latest/volcello_Ofx_'+model_name+'_historical_r1i1p1f1_gn.nc']
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list)

# EC-Earth3-CC
# -> OK!
#centre_name='EC-Earth-Consortium'
#model_name='EC-Earth3-CC'
#data_list=glob.glob(data_dir+'o2_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so_*_'+model_name+'_*')
#grid_list=glob.glob(data_dir+'volcello_*_'+model_name+'_*')
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list)

# UKESM
# -> OK!
#centre_name='MOHC'
#model_name='UKESM1-0-LL'
#data_dir='../../input_nc/UKESM1-0-LL/'
#data_list=glob.glob(data_dir+'o2_*_'+model_name+'_'+exp+'_*')
#t_list=glob.glob(data_dir+'thetao_*_'+model_name+'_'+exp+'_*')
#s_list=glob.glob(data_dir+'so_*_'+model_name+'_'+exp+'_*')
#grid_list=glob.glob('/badc/cmip6/data/CMIP6/CMIP/MOHC/UKESM1-0-LL/historical/r1i1p1f2/Omon/volcello/gn/latest/*')
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list)

# CMCC
# -> OK!
# -> OK!
#centre_name='CMCC'
#model_name='CMCC-ESM2'
#data_list=glob.glob(data_dir+data+'_*_'+model_name+'_*')
#t_list=glob.glob(data_dir+'thetao'+'_*_'+model_name+'_*')
#s_list=glob.glob(data_dir+'so'+'_*_'+model_name+'_*')
#grid_list=['/badc/cmip6/data/CMIP6/CMIP/'+centre_name+'/'+model_name+'/historical/r1i1p1f1/Ofx/volcello/gn/latest/volcello_Ofx_'+model_name+'_historical_r1i1p1f1_gn.nc']*2
#hfs.hackathon_wrapper(data_list,grid_list,data,exp,model_name,timeslices,output_dir,t_list,s_list)
