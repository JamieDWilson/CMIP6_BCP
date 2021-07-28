import xarray as xr
import numpy as np
import seawater as sw
import glob

def dim_info(data):
	   
# dim info()
# 
# Description: outputs dimension names of variables in file in consistent order. If gridname is not empty, also checks against grid variables 
# returning 'NaN' if different sizes.
#
# Inputs:
#	data - xarray array
#
# Output:
#	returned array - longitude name, latitude name, depth name 
#
# Author: Jamie D. Wilson
	
	# dimension info
	info=data.dims
    
	# longitude dimension name
	# -> add additional names if needed!
	if 'i' in data.dims:
		lon_var='i'
	elif 'lon' in data.dims:
		lon_var='lon'
	elif 'x' in data.dims:
		lon_var='x'
	elif 'nlon' in data.dims:
		lon_var='nlon'
	elif 'xh' in data.dims:
		lon_var='xh'
	else:
		lon_var='NaN'
		print(data.dims)
    
	# latitude dimension name
	# -> add additional names if needed!
	if 'j' in data.dims:
		lat_var='j'
	elif 'lat' in data.dims:
		lat_var='lat'
	elif 'y' in data.dims:
		lat_var='y'
	elif 'nlat' in data.dims:
		lat_var='nlat'
	elif 'yh' in data.dims:
		lat_var='yh'  
	else:
		lat_var='NaN'
		print(data.dims)

	# get depth dimensions if 3D
	# -> add additional names if needed!
	if len(info)>2:
		if 'depth' in data.dims:
			depth_var='depth'
		elif 'lev' in data.dims:
			depth_var='lev'
		elif 'olevel' in data.dims:
			depth_var='olevel'
		elif 'deptht' in data.dims:
			depth_var='deptht'
		else:
			depth_var='NaN'
			print(data.dims)
	else:
		depth_var='NaN'
            
	# output
	return [lon_var,lat_var,depth_var]
    
    ####### END ######


def compute_transfer_efficiency(pocdata,areadata,timeslice_times,model_name,output_dir):

# compute_transfer_efficiency()
#
# description: reads in expc and area files, extracts POC flux at two depths, calculates annual mean timeseries and timeslices of POC at two depths and ratio
#
# inputs:
#	pocdata - xarray dataset containing expc 
#	griddata - xarray dataset containing areacello
#	timeslice_times - pairs of years to average timeslice over, e.g., ['1960','1970'] or ['1960','1970','1980' '1990'] 
#	model_name - name of model for saving output
#	output_dir - relative path to output directory (include / at the end!)
#	flag_interpolate - True to linearly interpolate to specified depths or False to choose nearest depths on grid
#
# outputs:
#	exp_timeseries_* - annual total POC flux at upper depth (Pg C year-1)
#	deep_timeseries_* - annual total POC flux at lower depth (Pg C year-1)
#	te_timeseries_* - annual area-weighted mean ratio of flux at depths (unitless)
#	exp_[year1]_[year2]_* - POC flux at upper depth averaged over timeslice years (mol m-2 yr-1)
#|      deep_[year1]_[year2]_* - POC flux at lower depth averaged over timeslice years (mol m-2 yr-1)
#|      te_[year1]_[year2]_* - ratio of POC fluxes at depths averaged over timeslice years (unitless)
#
# author: Jamie D. Wilson (Francisco de Melo Virissimo; Rui Ying; Markus Adloff)
# 
# To do:
#	1) check interpolation works - DONE 
#	2) automatic dimensions / change to standard dimension - NEEDS TESTING
#	3) expc variable name disappearing in timeseries - DONE
#	4) check weighting code is correct - DONE - maybe local persmissions issue?
#	5) netcdf overwrite permissions - DONE

	######################################
        ###### hard coded user options #######
	flag_interpolate=True
	exp_horizon=100 # export horizon (m)
	poc_horizon=1000 # te depth horizon (m)
	sec_to_yr=60.0*60.0*24.0*365.0 # seconds to year
	mol_to_g_C=12.01 # mol to gram for carbon
	g_to_Pg=1.0/1.0e15 # grams to Petagrams
	#####################################

	poc=pocdata.expc
	area=areadata.areacello

	# get info to force consistent renaming of dimensions/coordinates
	info=dim_info(poc)
	poc=poc.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})
	info=dim_info(area)
	area=area.rename({info[0]:'i',info[1]:'j'})
        
	# get arrays
	if flag_interpolate:

		# depth levels to interpolate between
		exp_lev_flag=False
		exp_lev_1=poc.sel(lev=exp_horizon,method='nearest').lev.data
		exp_lev_2=poc.sel(lev=exp_horizon,method='backfill').lev.data
		if exp_lev_2 == exp_lev_1:
			exp_lev_2=poc.sel(lev=exp_horizon,method='pad').lev.data
			exp_lev_flag=True
	
		poc_lev_flag=False
		poc_lev_1=poc.sel(lev=poc_horizon,method='nearest').lev.data
		poc_lev_2=poc.sel(lev=poc_horizon,method='backfill').lev.data
		if poc_lev_2 == poc_lev_1:
			poc_lev_2 = poc.sel(lev=poc_horizon,method='pad').lev.data
			poc_lev_flag=True

		exp_levels=[exp_lev_1,exp_lev_2]
		if exp_lev_1==exp_horizon and exp_lev_2==exp_horizon:
			 exp_scalar=0.0
		else:
			exp_scalar=(exp_horizon-min(exp_levels)) / (max(exp_levels)-min(exp_levels))

		poc_levels=[poc_lev_1,poc_lev_2]
		if poc_lev_1==poc_horizon and poc_lev_2==poc_horizon:
			poc_scalar=0.0
		else:
			poc_scalar=(poc_horizon-min(poc_levels)) / (max(poc_levels)-min(poc_levels))

		# get data at levels - upper
		print('>>> interpolating poc data between',min(exp_levels),'and',max(exp_levels),'to',exp_horizon,'m')
		flux_exp_1 = poc.sel(lev=exp_horizon, method="nearest").compute() # export flux at nearest level above z = depth
		if exp_lev_flag:
			flux_exp_2 = poc.sel(lev=exp_horizon, method="pad").compute() # export flux at other level
		else:
			flux_exp_2 = poc.sel(lev=exp_horizon, method="backfill").compute() 
		
		if exp_levels[0] < exp_levels[1]:
			exp = (flux_exp_1 * (1.0 - exp_scalar)) + (flux_exp_2 * exp_scalar)
		else:
			exp = (flux_exp_1 * exp_scalar) + (flux_exp_2 *(1.0 - exp_scalar))
		
                # get data at levels - lower
		print('>>> interpolating poc data between',min(poc_levels),'and',max(poc_levels),'to',poc_horizon,'m')
		flux_poc_1 = poc.sel(lev=poc_horizon, method="nearest").compute() # export flux at nearest level above z = depth
		if poc_lev_flag:
                        flux_poc_2 = poc.sel(lev=poc_horizon, method="pad").compute() # export flux at other level
		else:
                        flux_poc_2 = poc.sel(lev=poc_horizon, method="backfill").compute()
                
		if poc_levels[0] < poc_levels[1]:
			deep = (flux_poc_1 * (1.0 - poc_scalar)) + (flux_poc_2 * poc_scalar)
		else:
			deep = (flux_poc_1 * poc_scalar) + (flux_poc_2 *(1.0 - poc_scalar))

	else:
		exp=poc.sel(lev=exp_horizon,method='nearest')
		deep=poc.sel(lev=poc_horizon,method='nearest')
		print('>>> selecting nearest depth levels:',poc.sel(lev=exp_horizon,method='nearest').lev.data,poc.sel(lev=poc_horizon,method='nearest').lev.data)
	
	print('>>> calculating annual averages')
        ### annual averaging
        # compue here to save memory in following steps
	exp=exp.groupby('time.year').mean('time').compute()
	deep=deep.groupby('time.year').mean('time').compute()

	print('>>> calculating transfer efficiency')
        ### calculate transfer efficiency
	te=deep/exp

        # calculate e-folding depth
        # TO DO? 

	print('>>> output transfer efficiency timeseries')
        # convert to global mean timeseries
	# weight by area
	weights = area/area.sum(dim=['i','j'])
	weights = weights.fillna(0)
	weighted_te = te.weighted(weights)
               
	# convert to global mean timeseries	
	#te_timeseries=weighted_te.sum(dim=['i','j'])
	te_timeseries=weighted_te.mean(dim=['i','j'])
	# output to netCDF
	te_timeseries=te_timeseries.rename('expc')
	#te_timeseries.to_netcdf(output_dir+'te_timeseries_'+expc_name,mode='w')
	
        # get timeslice data
	years=te.year
	for n in range(0,len(timeslice_times),2):

		if (years==int(timeslice_times[n])).any() and (years==int(timeslice_times[n+1])).any():
			print('>>> output transfer efficiency timeslices:',int(timeslice_times[n]),'and',int(timeslice_times[n+1]))
			te_timeslice=te.sel(year=slice(timeslice_times[n], timeslice_times[n+1])).mean(dim='year')
			te_timeslice.to_netcdf(output_dir+'te_'+timeslice_times[n]+'_'+timeslice_times[n+1]+'_'+model_name+'.nc',mode='w')
	
	print('>>> output timseries of export and poc at depth')
	# output the timeseries of total POC export at 100m
	exp_timeseries = ((exp*area)*sec_to_yr*mol_to_g_C*g_to_Pg).sum(dim=['i','j'])
	exp_timeseries=exp_timeseries.rename('exp')
	#exp_timeseries.to_netcdf(output_dir+'exp_timeseries_'+expc_name,mode='w')

	deep_timeseries = ((deep*area)*sec_to_yr*mol_to_g_C*g_to_Pg).sum(dim=['i','j'])
	deep_timeseries=deep_timeseries.rename('deep')
	#deep_timeseries.to_netcdf(output_dir+'deep_timeseries_'+expc_name,mode='w')

        # get export timeslice data
	years=exp.year
	total_exp=np.zeros(int(len(timeslice_times)/2))
	count=0
	for n in range(0,len(timeslice_times),2):

		if (years==int(timeslice_times[n])).any() and (years==int(timeslice_times[n+1])).any():
			print('>>> output timeslices of export and poc at depth:',int(timeslice_times[n]),'and',int(timeslice_times[n+1]))
			exp_timeslice=exp.sel(year=slice(timeslice_times[n],timeslice_times[n+1])).mean(dim='year')
			exp_timeslice.to_netcdf(output_dir+'exp_'+timeslice_times[n]+'_'+timeslice_times[n+1]+'_'+model_name+'.nc',mode='w')
			
			deep_timeslice=deep.sel(year=slice(timeslice_times[n],timeslice_times[n+1])).mean(dim='year')
			deep_timeslice.to_netcdf(output_dir+'deep_'+timeslice_times[n]+'_'+timeslice_times[n+1]+'_'+model_name+'.nc',mode='w')

			total_exp[count] = ((exp_timeslice*area)*sec_to_yr*mol_to_g_C*g_to_Pg).sum(dim=['i','j']).compute().data
			count=count+1
	
	# return total export from timeslices if collecting in looping script
	#return [te_timeseries,exp_timeseries,deep_timeseries,total_exp]
	
	timeseries=xr.merge([te_timeseries,exp_timeseries,deep_timeseries],compat='override')
	return [timeseries,total_exp]



def compute_export_production(pocdata,areadata,timeslice_times,model_name,output_dir):

# compute_export_production()
#
# description: reads in epc100 and area files, calculates annual total timeseries and timeslices of POC export
#
# inputs:
#|      pocdata - xarray dataset containing epc100
#|      griddata - xarray dataset containing areacello
#|      timeslice_times - pairs of years to average timeslice over, e.g., ['1960','1970'] or ['1960','1970','1980' '1990'] 
#|      model_name - name of model when saving output
#	output_dir - relative path to output directory (include / at the end!)
#
# outputs:
#|      epc100_timeseries_* - annual total POC export (Pg C year-1)
#|      epc100_[year1]_[year2]_* - POC averaged over timeslice years (mol m-2 yr-1)
#
# author: Jamie D. Wilson (Katie Sieradzan, Fraser Goldsworth)

	######################################
        ###### hard coded user options #######
	sec_to_yr=60.0*60.0*24.0*365.0 # seconds to year
	mol_to_g_C=12.01 # mol to gram for carbon
	g_to_Pg=1.0/1.0e15 # grams to Petagrams
	####################################

	# testing generic wrapper
	poc=pocdata.epc100
	area=areadata.areacello
	
	# get info to force consistent renaming of dimensions/coordinates
	info=dim_info(poc)
	poc=poc.rename({info[0]:'i',info[1]:'j'})
	area=area.rename({info[0]:'i',info[1]:'j'})
	
	print('>>> calculating annual averages')
        ### annual averaging
        # compue here to save memory in following steps
	poc=poc.groupby('time.year').mean('time').compute()
	
	print('>>> output timeseries')
	# convert to global total timeseries    
	poc_timeseries=((poc*area)*sec_to_yr*mol_to_g_C*g_to_Pg).sum(dim=['i','j'])
	# output to netCDF
	poc_timeseries=poc_timeseries.rename('epc100')
	poc_timeseries=poc_timeseries.to_dataset()
	#poc_timeseries.to_netcdf(output_dir+'epc100_timeseries_'+model_name+'.nc',mode='w')

	# get timeslice data
	years=poc.year
	total_poc=np.zeros(int(len(timeslice_times)/2))
	count=0
	for n in range(0,len(timeslice_times),2):

		if (years==int(timeslice_times[n])).any() and (years==int(timeslice_times[n+1])).any():
			print('>>> output selected timeslices:',int(timeslice_times[n]),'and',int(timeslice_times[n+1]))
			poc_timeslice=poc.sel(year=slice(timeslice_times[n], timeslice_times[n+1])).mean(dim='year')
			poc_timeslice.to_netcdf(output_dir+'epc100_'+timeslice_times[n]+'_'+timeslice_times[n+1]+'_'+model_name+'.nc',mode='w')
			total_poc[count] = ((poc_timeslice*area)*sec_to_yr*mol_to_g_C*g_to_Pg).sum(dim=['i','j']).compute().data
			count=count+1

	# return total export from timeslices if collecting in looping script
	return [poc_timeseries,total_poc]
        
###### END ##### 


def compute_AOU_Csoft(dataset,griddata,timeslice_times,model_name,output_dir):
# compute_AOU_Csoft()
#
# description: takes in dataset with o2, temperature, salinity and volume, calculates O2 saturation, calculates annual mean/total timeseries and timeslices of AOU & Csoft
#
# inputs:
#|      daaset - xarray datset containing o2, thetao, so, (and volume if 4D)
#|      griddata - xarray dataset containing volume
#|      timeslice_times - pairs of years to average timeslice over, e.g., ['1960','1970'] or ['1960','1970','1980' '1990'] 
#|      model_name - name of model when saving output
#	output_dir - relative path to output directory (include / at the end!)
#
# outputs:
#|      AOU_[year1]_[year2]_* - AOU averaged over timeslice years (mol m-3)
#|      Csoft_[year1]_[year2]_* - Csoft averaged over timeslice years (Pg C)
#|      returned array - Csoft_timeseries and total_Csoft (Pg C)
#
# author: Jamie D. Wilson (adapted from code by Anna Katavouta, Emily Vosper)


	######################################
	###### hard coded user options #######
	sec_to_yr=60.0*60.0*24.0*365.0 # seconds to year
	mol_to_g_C=12.01 # mol to gram for carbon
	g_to_Pg=1.0/1.0e15 # grams to Petagrams
	redfield_C_O=106/170
	####################################

	o2=dataset.o2
	t=dataset.thetao
	s=dataset.so
	if len(griddata.volcello.dims)>3:
		volume=dataset.volcello
	else:
		volume=griddata.volcello

	# get info to force consistent renaming of dimensions/coordinates
	info=dim_info(o2)
	o2=o2.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})

	info=dim_info(volume)
	volume=volume.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})

	info=dim_info(t)
	t=t.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})

	info=dim_info(s)
	s=s.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})

	print('>>> calculating annual averages...')
	### annual averaging
	# compue here to save memory in following steps
	print('>> ...for O2')
	o2=o2.groupby('time.year').mean('time').compute()
	print('>> ...for T')
	t=t.groupby('time.year').mean('time').compute()
	print('>> ...for S')
	s=s.groupby('time.year').mean('time').compute()
	if len(volume.dims)>3:
		print('>> ...for volume (because it is varying with time)')
		volume=volume.groupby('time.year').mean('time').compute()

	print('>>> calculate oxygen saturation and AOU')
	# oxygen at saturation 
	oxygen_sat = compute_o2sat(t,s)
	dens = sw.eos80.dens0(s,t)
	# conver from to ml l-1 to mol m-3
	oxygen_sat = oxygen_sat * dens
	#oxygen_sat_per =  (oxygen/oxygen_sat)*100.0
	AOU = (oxygen_sat - o2).compute()

	print('>>> calculate total Csoft')
	# volume averaging nd calculate global total time-series
	Csoft=((AOU*volume) * mol_to_g_C * g_to_Pg *redfield_C_O).compute()
	
	print('>>> output timeseries of Csoft')
	Csoft_timeseries = Csoft.sum(dim=['lev','i','j'])
	Csoft_timeseries=Csoft_timeseries.rename('o2')
	Csoft_timeseries=Csoft_timeseries.to_dataset()
	#Csoft_timeseries.to_netcdf(output_dir+'Csoft_timeseries_'+o2_name,mode='w')

	years=Csoft.year
	total_Csoft=np.zeros(int(len(timeslice_times)/2))
	count=0
	for n in range(0,len(timeslice_times),2):
		if (years==int(timeslice_times[n])).any() and (years==int(timeslice_times[n+1])).any():
			print('>>> output selected timeslices of Csoft and AOU:',int(timeslice_times[n]),'and',int(timeslice_times[n+1]))
			Csoft_timeslice=Csoft.sel(year=slice(timeslice_times[n], timeslice_times[n+1])).mean(dim='year')
			Csoft_timeslice=Csoft_timeslice.rename('o2')
			Csoft_timeslice.to_netcdf(output_dir+'Csoft_'+timeslice_times[n]+'_'+timeslice_times[n+1]+'_'+model_name+'.nc',mode='w')
			total_Csoft[count] = Csoft_timeslice.sum(dim=['lev','i','j']).compute().data
			AOU_timeslice=AOU.sel(year=slice(timeslice_times[n], timeslice_times[n+1])).mean(dim='year')
			AOU_timeslice=AOU_timeslice.rename('o2')
			AOU_timeslice.to_netcdf(output_dir+'AOU_'+timeslice_times[n]+'_'+timeslice_times[n+1]+'_'+model_name+'.nc',mode='w')
			count=count+1
	

	# return total export from timeslices if collecting in looping script
	return [Csoft_timeseries,total_Csoft]

###### END #####



def compute_o2sat(t,s):

# compute_o2sat()
# Description: calculates oxygen saturation from temperature and salinity .Oxygen saturation value is the volume of oxygen gas absorbed from humidity-saturated
# air at a total pressure of one atmosphere, per unit volume of the liquid at the temperature of measurement
#
# Inputs: 
#	t - temperature (units?)
#	s - salinity (units?)
#
# Outputs:
#	oxsat - oxygen saturation (ml l-1)
#
# Author: Anna Katavouta, adapted by Jamie D. Wilson

	oA0=  2.00907;
	oA1=  3.22014;
	oA2=  4.05010;
	oA3=  4.94457;
	oA4= -2.56847E-1;
	oA5=  3.88767;
	oB0= -6.24523E-3;
	oB1= -7.37614E-3;
	oB2= -1.03410E-2;
	oB3= -8.17083E-3;
	oC0= -4.88682E-7;

	aTT = 298.15-t;
	aTK = 273.15+t;
	aTS = np.log(aTT/aTK);
	aTS2= aTS*aTS;
	aTS3= aTS2*aTS;
	aTS4= aTS3*aTS;
	aTS5= aTS4*aTS;

	ocnew= np.exp(oA0 + oA1*aTS + oA2*aTS2 + oA3*aTS3 + oA4*aTS4 + oA5*aTS5 + \
                  s*(oB0 + oB1*aTS + oB2*aTS2 + oB3*aTS3) + oC0*(s*s));

	# to units mol/kg
	oxsat = ocnew/22391.6;
	return oxsat
 
####### END #########


def preprocess(dir_list,timeslice_years):

	# preprocess()
	# 
	# Description: loads in years of data in dir_list, checks where timelice_years fall within, sets
	# a flag variable used to determine whether to run analysis function or concanenate output
	#
	# Inputs:
	# dir_list - list of file paths to CMIP6 data
	# timeslice_years - list of year-pairs that define timeslice for averaging, e.g., ['1960','1970']
	#
	# Author: Jamie D. Wilson

	# list output files for specific model in specific directory
	n_files=len(dir_list)
	
	if timeslice_years:

		# preprocess runs - concatenate runs so that timeslices are fully within for averaging
		# 1/2) find start/end years of each run, check if timeslice years exist within range 
		startend=np.zeros((len(dir_list),2))

		for n in range(0,len(dir_list)):
			data=xr.open_mfdataset(dir_list[n])
			year_end=max(data.time.dt.year.data)
			year_start=min(data.time.dt.year.data)
			for nn in range(0,len(timeslice_years),2):
				if int(timeslice_years[nn])>=year_start and int(timeslice_years[nn])<=year_end:
					startend[n,0]=1
				if int(timeslice_years[nn+1])>=year_start and int(timeslice_years[nn+1])<=year_end:
					startend[n,1]=1
					startend[n,0]=0

		# 2/2) generate flags for compute or concatenate 
		compute_flag=[]
		n=0
		main_flag=True
		while main_flag:
			if startend[n,0]==1. and startend[n,1]==0.:
				compute_flag.append(False)
				inner_flag=True
				nn=n+1 # initiate next counter
				while inner_flag:
					if startend[nn,0]==0. and startend[nn,1]==1.:
						compute_flag.append(True)
						n=nn+1 # reset original counter
						inner_flag=False # quit this while loop
					else:
						compute_flag.append(False)
						nn=nn+1 # continue
			else:
				compute_flag.append(True)
				n=n+1

			# stop while loop if exceeding number of files
			if n>n_files-1:
				main_flag=False
	
	else: 
		compute_flag=[True]*n_files

	return compute_flag


def hackathon_wrapper(data_list,grid_list,data_name,exp_name,model_name,timeslice_years,output_dir,t_list=[],s_list=[],thk_list=[]):

	# hackathon_wrapper()
	#
	# Description: Finds any CMIP6 experiments that need concatenating for averaging, runs analysis function
	# and concatenates timeseries output into single file
	#
	# Inputs:
	# data_list - list of paths to CMIP6 output
	# grid_list - list of paths to associated grid output (area or volume)
	# exp_name - name of specific experiment, e.g., 'historical', used to name output files
	# model_name - name of model, e.g., 'UKESM', used to name output files
	# timeslice_years - list of year-pairs to average timeslices over, e.g., ['1960','1970','2000','2010']
	# output_dir - output directory path (end with '/')
	# t_list - optional list ofinput for temperature files
	# s_list - optional list of paths for salinty files
	#
	# Author: Jamie D. Wilson


	# check inputs
	if not data_list:
		raise NameError('data_list is empty')
	if not grid_list:
		raise NameError('grid_list is empty')

	# assign specific compute function to common function name
	if data_name=='epc100':
		compute=compute_export_production
		grid_name='areacello'
	elif data_name=='expc':
		compute=compute_transfer_efficiency
		grid_name='areacello'
	elif data_name=='o2':
		compute=compute_AOU_Csoft
		grid_name='volcello'
	else: 
		print('unknown data name')

	# preprocessing - figure out if inputs need concatenating for averaging
	compute_flag=preprocess(data_list,timeslice_years)

	# initialise main timeseries output file
	main_timeseries=xr.Dataset(dict(data=(['year'],[])),coords=dict(year=[]))
	main_timeseries=main_timeseries.rename({'data':data_name}) 
	if data_name=='expc':
		# add additional timseries outputs
		main_timeseries['exp']=xr.DataArray(dims=['year'],coords=dict(year=[])).rename('exp')
		main_timeseries['deep']=xr.DataArray(dims=['year'],coords=dict(year=[])).rename('deep')

	# load in grid variable - if 4D then this is not used 
	grid=xr.open_mfdataset(grid_list[0])

	# convert area to volume if inputs are o2 and area
	if thk_list:
		thk=xr.open_mfdataset(thk_list[0])
		grid=area_to_volume(grid,thk)

	# load in first set of output(s) 
	#if t_list or s_list:
	#	# AOU function takes more variables
	#	if len(grid[grid_name].dims)>3:
	#		data_concat=xr.open_mfdataset([data_list[0],t_list[0],s_list[0],grid_list[0]])
	#		data_concat=check_inputs(data_concat,data_name)
	#	else:
	#		data_concat=xr.open_mfdataset([data_list[0],t_list[0],s_list[0]])
	#		data_concat=check_inputs(data_concat,data_name)	
	#else:
	#	data_concat=xr.open_mfdataset(data_list[0])
	#	data_concat=check_inputs(data_concat,data_name)
	data_concat=load_data(data_name,grid_name,grid,data_list,t_list,s_list,grid_list,thk_list,0)

	for n in range(1,len(data_list)):
	
		# load in second set of output(s)
	#	if t_list or s_list:
	#		if len(grid[grid_name].dims)>3:
	#			data_2=xr.open_mfdataset([data_list[n],t_list[n],s_list[n],grid_list[n]])
	#			data_2=check_inputs(data_2,data_name)
	#		else:
	#			data_2=xr.open_mfdataset([data_list[n],t_list[n],s_list[n]])
	#			data_2=check_inputs(data_2,data_name)
	#	else:
	#		data_2=xr.open_mfdataset(data_list[n])
	#		data_2=check_inputs(data_2,data_name)
		data_2=load_data(data_name,grid_name,grid,data_list,t_list,s_list,grid_list,thk_list,n)

		if compute_flag[n-1]:
			print('computing file:',data_list[n-1],'for',data_name)
			# compute step
			returned_dataset=compute(data_concat,grid,timeslice_years,model_name,output_dir)
			single_timeseries=returned_dataset[0]
			#print(data_processed.dims,data_out.dims)
			# post-processing - combine all timeseries outputs
			main_timeseries=xr.concat([main_timeseries,single_timeseries],dim='year')
			# reset 
			data_concat=data_2
		else:
			# concatenate files together
			print('concatenating file:',data_list[n])
			data_concat=xr.concat([data_concat,data_2],dim='time')

	# compute the remaining data
	print('computing file:',data_list[-1],' for',data_name)
	returned_dataset=compute(data_concat,grid,timeslice_years,model_name,output_dir)
	single_timeseries=returned_dataset[0]
	main_timeseries=xr.concat([main_timeseries,single_timeseries],dim='year')

	# write final timeseries file
	print('writing final output to:',output_dir+data_name+'_'+exp_name+'_'+model_name+'.nc')
	main_timeseries.to_netcdf(output_dir+data_name+'_'+exp_name+'_'+model_name+'.nc',mode='w')


def check_inputs(dataset,dataname):

	# get dimension info
	info=dim_info(dataset[dataname])

	if len(info)>2:
		depth_dim_name=info[2]
		
		if dataset[depth_dim_name].units=='centimeters':
			dataset[depth_dim_name]=dataset[depth_dim_name]/100

	return dataset


def area_to_volume(area,thickness):

	print('>>> converting areacello to volcello')
	# correct cm to m in depth
	thickness=check_inputs(thickness,'thkcello')
	
	# calculate volume
	volume=area.areacello*thickness.thkcello

	# output to dataset with correct variable name
	volume=volume.rename('volcello')
	volume=volume.to_dataset()

	return volume


def load_data(data_name,grid_name,grid,data_list,t_list,s_list,grid_list,thk_list,n):

	# load in output(s) 
	if t_list or s_list:
		
		# AOU function takes more variables
		data_out=xr.open_mfdataset([data_list[n],t_list[n],s_list[n]])
		
		# add volume if varying with time
		if len(grid[grid_name].dims)>3:
			
			# convert from area to volume if necessary
			if thk_list:
				tmp=xr.open_mfdataset(grid_list[0]) # single area file
				thk=xr.open_mfdataset(thk_list[n])
				tmp=area_to_volume(tmp,thk)
			else:
				tmp=xr.open_mfdataset(grid_list[n]) # multiple volume files

			# merge with data_concat
			data_out=data_out.merge(tmp)
			
	else:
		data_out=xr.open_mfdataset(data_list[n])
	
	data_out=check_inputs(data_out,data_name)

	return data_out
