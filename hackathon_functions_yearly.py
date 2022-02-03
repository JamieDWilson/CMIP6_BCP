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

def compute_AOU_Csoft(dataset_o2,dataset_T,dataset_S,dataset_V,model_name,nstart,nend,output_dir):
# compute_AOU_Csoft()
#
# description: takes in dataset with o2, temperature, salinity and volume, calculates O2 saturation, calculates annual mean/total timeseries and timeslices of AOU & Csoft
#
# inputs:
#|      dataset - xarray datset containing o2, thetao, so, (and volume if 4D)
#|      model_name - name of model when saving output
#	output_dir - relative path to output directory (include / at the end!)
#
# outputs:
#|      returned array - Csoft_timeseries (monthly or yearly)
#
# author: Jamie D. Wilson (adapted from code by Anna Katavouta, Emily Vosper)


	######################################
	###### hard coded user options #######
	sec_to_yr=60.0*60.0*24.0*365.0 # seconds to year
	mol_to_g_C=12.01 # mol to gram for carbon
	g_to_Pg=1.0/1.0e15 # grams to Petagrams
	redfield_C_O=106/170
	####################################

	o2=dataset_o2.o2[nstart:nend,:,:,:]
	t=dataset_T.thetao[nstart:nend,:,:,:]
	s=dataset_S.so[nstart:nend,:,:,:]
	if len(dataset_V.volcello.dims)>3:
		volume=dataset_V.volcello[nstart:nend,:,:,:]
	else:
		volume=dataset_V.volcello

	# get info to force consistent renaming of dimensions/coordinates
	info=dim_info(o2)
	o2=o2.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})

	info=dim_info(volume)
	volume=volume.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})

	info=dim_info(t)
	t=t.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})

	info=dim_info(s)
	s=s.rename({info[0]:'i',info[1]:'j',info[2]:'lev'})
        
        # If you want annual averages uncomments these lines
	print('>>> calculating annual averages...')
	### annual averaging
	## compue here to save memory in following steps
	print('>> ...for O2')
	o2=o2.groupby('time.year').mean('time').compute()
	#print('>> ...for T')
	t=t.groupby('time.year').mean('time').compute()
	#print('>> ...for S')
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
	Csoft_timeseries=Csoft_timeseries.rename('Csoft')
	Csoft_timeseries=Csoft_timeseries.to_dataset()
	#Csoft_control = xr.Dataset(
#		{
#			"Csoft_control": (["time"], Csoft_timeseries.Csoft),
#			"time_years": (["time"], time_years),
#		},
#		coords={
#			"time" : time_all[:],
#		},
#	)
#
#	Csoft_control.to_netcdf(output_dir+'Csoft_'+exp+'_'+model_name,mode='w')
	#Csoft_timeseries.to_netcdf(output_dir+'Csoft_timeseries_'+o2_name,mode='w')
        
	return Csoft_timeseries
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


def compute_transfer_efficiency(pocdata,areadata,timeslice_times,datch,model_name,exp_name,nstart,nend,output_dir):

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
#       datch - for modifying year in control simulation to match actual years like 2050 etc.
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

	poc=pocdata.expc[nstart:nend,:,:,:]
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
	timeslices_write=[element+datch for element in timeslice_times]
	years=te.year
	for n in range(0,len(timeslice_times),2):
		if (years==int(timeslice_times[n])).any() and (years==int(timeslice_times[n+1])).any():
			print('>>> output transfer efficiency timeslices:',int(timeslices_write[n]),'and',int(timeslices_write[n+1]))
			te_timeslice=te.sel(year=slice(str(timeslice_times[n]), str(timeslice_times[n+1]))).mean(dim='year')
			te_timeslice.to_netcdf(output_dir+'te_'+str(timeslices_write[n])+'_'+str(timeslices_write[n+1])+'_'+model_name+'_'+exp_name+'.nc',mode='w')
	
	print('>>> output timseries of export and poc at depth')
	# output the timeseries of total POC export at 100m
	exp_timeseries = ((exp*area)*sec_to_yr*mol_to_g_C*g_to_Pg).sum(dim=['i','j'])
	exp_timeseries=exp_timeseries.rename('exp')
	#exp_timeseries.to_netcdf(output_dir+'exp_timeseries_'+expc_name/expc_name,mode='w')

	deep_timeseries = ((deep*area)*sec_to_yr*mol_to_g_C*g_to_Pg).sum(dim=['i','j'])
	deep_timeseries=deep_timeseries.rename('deep')
	#deep_timeseries.to_netcdf(output_dir+'deep_timeseries_'+expc_name,mode='w')

        # get export timeslice data
	years=exp.year
	total_exp=np.zeros(int(len(timeslice_times)/2))
	count=0
	for n in range(0,len(timeslice_times),2):

		if (years==int(timeslice_times[n])).any() and (years==int(timeslice_times[n+1])).any():
			print('>>> output timeslices of export and poc at depth:',int(timeslices_write[n]),'and',int(timeslices_write[n+1]))
			exp_timeslice=exp.sel(year=slice(str(timeslice_times[n]),str(timeslice_times[n+1]))).mean(dim='year')
			exp_timeslice.to_netcdf(output_dir+'exp_'+str(timeslices_write[n])+'_'+str(timeslices_write[n+1])+'_'+model_name+'_'+exp_name+'.nc',mode='w')
			
			deep_timeslice=deep.sel(year=slice(str(timeslice_times[n]),str(timeslice_times[n+1]))).mean(dim='year')
			deep_timeslice.to_netcdf(output_dir+'deep_'+str(timeslices_write[n])+'_'+str(timeslices_write[n+1])+'_'+model_name+'_'+exp_name+'.nc',mode='w')

			total_exp[count] = ((exp_timeslice*area)*sec_to_yr*mol_to_g_C*g_to_Pg).sum(dim=['i','j']).compute().data
			count=count+1
	
	# return total export from timeslices if collecting in looping script
	#return [te_timeseries,exp_timeseries,deep_timeseries,total_exp]
	
	timeseries=xr.merge([te_timeseries,exp_timeseries,deep_timeseries],compat='override')
	return [timeseries,total_exp]

