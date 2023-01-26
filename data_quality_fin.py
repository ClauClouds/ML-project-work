#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import datetime as dt
from datetime import datetime, timedelta
from suntime import Sun, SunTimeException
from datetime import datetime, timezone
import numpy as np
import glob
import os
import xarray as xr
from pathlib import Path
import gzip
import shutil


# In[ ]:


####################### all functions ###########################

def sunlight_juelich(day, ds, time2):
    '''
    For the day it selects the timesteps with daylight (one hour after sunrise and one hour after sunrise)
    !! time dimension from dataset needs to be called time or time10s and time60s!!
    
    time2: use for datasets from the solarpanel, for the rest put time2 = None'''

    #coordinates of Juelich
    latitude = 50.9224
    longitude = 6.3639
    sun = Sun(latitude, longitude)
    
    #getting the time of sunrise +1 hour and the time of sunset - 1 hour. (result has information of timezone)
    sur_tzd = sun.get_local_sunrise_time(day) + timedelta(hours=1)
    sus_tzd = sun.get_local_sunset_time(day) - timedelta(hours=1)

    #delete timezonedesignator (Data from sup_joys are in UTC)
    sur = sur_tzd.replace(tzinfo=None)
    sus = sus_tzd.replace(tzinfo=None)
    if time2 == None:
        ds_light = ds.sel(time=slice(sur, sus))
    else:
        ds_light = ds.sel(time10s=slice(sur, sus), time60s=slice(sur, sus))
   
    return ds_light



def merge_datasets(startdate, enddate, outpath, outname):
    '''
    What does it?
    Step 4 from the Instructions:
    
    4. Saving data in array dataset instead of pandas data frame
        1. Iterate on a list of files from input
        2. Read also SW and LW downwelling, (as well as direct/diffuse SW radiation)
    
    
    -> In a timeperiod (between startdate and enddate) opens the daily datasets, puts all data (LW, LW_error, SW, SW_error) 
    into one dataset and selects the timstamps with sunlight.
    Saves it as a dataset in a given name and path. It creates new folders for every year if they do not already exist.
    
    
    input:
        - startdate: <string>: e.g.: '2020-01-01' 
        - enddate:   <string>: e.g.: '2020-03-02'
        - outpath:   <string>: where to save the nc-file eg.: 'full_datasets'
        - outname:   <string>: e.g.: 'all_data'  automaticlly appends the date of the day and an .nc at the end
        
        
    output:
        - saves a dataset (rlds, rlds_error, rsds, rsds_error), with the timesteps of the sunlight    
    ''' 
    #Range of the datas
    all_dates = pd.date_range(start=startdate,end=enddate)
    

    for day in all_dates:
        
        #get the day to string, to open the right file
        date = day.strftime('%Y%m%d') + '000000.nc'
        y = day.strftime('%Y')

        #Open Files for LW and shortwave data and put them in one dataset
         #getting the names from the data. Some names have ..._v00_... and some ..._v01...
        fileName_lw = glob.glob('/data/hatpro/jue/hdcp2/radiation_hdcp2/'+ y +'/sups_joy_pyrg00_l1_rlds_v*_'+ date )
        #if there is no file for the day it goes to the next day
        if not fileName_lw:
            print(date + ' no LW-file for this day found')
            continue

        #open files 
        fileObj_lw = Path(str(fileName_lw[0]))
        if fileObj_lw.is_file() == True:
            ds_lw = xr.open_dataset(str(fileName_lw[0]),
                                   drop_variables = ['lon', 'lat'])

        else: 
            print('File not found', fileName_lw[0])
            continue

        #getting the names from the data. Some names have ..._v00_... and some ..._v01...
        fileName_sw = glob.glob('/data/hatpro/jue/hdcp2/radiation_hdcp2/'+ y +'/sups_joy_pyr00_l1_rsds_v*_'+ date )
        #if there is no file for the day it goes to the next day
        if not fileName_sw:
            print(date + ' no SW-file for this day found')
            continue

        #open files 
        fileObj_sw = Path(str(fileName_sw[0]))
        if fileObj_sw.is_file() == True:
            ds_sw = xr.open_dataset(str(fileName_sw[0]))
        else: 
            print('File not found', fileName_sw[0])
            continue

        #put the sw data to the lw data
        ds_lw['rsds'] = ds_sw['rsds']
        ds_lw['rsds_error'] = ds_sw['rsds_error']
        ds = ds_lw
        
        #data with daylight
        ds_light = sunlight_juelich(day, ds, None)

        #checking if path/folder already exists
        MYDIR = (outpath + '/' + y)
        CHECK_FOLDER = os.path.isdir(MYDIR)

        # If folder doesn't exist, then create it.
        if not CHECK_FOLDER:
            os.makedirs(MYDIR)
            print("created folder : ", MYDIR)

        #saving it as a new dataset 
        save = outpath +'/'+ y + '/'+ outname + date 
        ds_light.to_netcdf(save)  





    


def data_quality(startdate, enddate, inpath, inname, outpath, outname):
    '''
    - you can choose a range of time, where it selects the datsets of the days where BOTH (LW and SW) Variables exists.
        (it uses the datasets which were created by the funct: merge_datasets)
    - checks the datasets for nan values (deletes the timesteps), unrealistic values (10<lw<700, -5<sw1500), and for gaps (non existent timesteps)
    - saves the gaps in a panda dataframes (starttime, endtime, length, how many timesteps), name of df e.g.: gaps_20180504.csv (just if there are gaps for that day)
    
     input:
        - startdate: <string>: e.g.: '2020-01-01' 
        - enddate:   <string>: e.g.: '2020-03-02'
        - inname:    <string>: e.g.: 'all_data'  automaticlly appends the date of the day and an .nc at the end
        - inpath:    <string>: where to save the nc-file eg.: 'full_datasets'
        - outpath:   <string>: where to save the csv-file eg.: 'data_quality'
        - outname:   <string>: name of the csv file (date is added): eg.: 'gaps' --> gaps_20200101
        
     output: 
        - saves gaps in a dataframe (starttime (the first missing), endtime (the last missing), length (counting one missing step as 5 sec), how many timesteps) as csv files (for each day)
    '''
    #Range of the datas
    all_dates = pd.date_range(start=startdate,end=enddate)

    #coordinates of Juelich
    latitude = 50.908546
    longitude = 6.413536
    sun = Sun(latitude, longitude)

    for day in all_dates:

        #get the day to string, to open the right file
        date = day.strftime('%Y%m%d') + '000000.nc'
        d = day.strftime('%d')
        m = day.strftime('%m')
        y = day.strftime('%Y')


        #Open Files for LW and shortwave data and put them in one dataset
         #getting the names from the data. Some names have ..._v00_... and some ..._v01...
        #fileName ='full_datasets/' + y + '/all_data' + date
        fileName = inpath +'/'+ y + '/'+ inname + date 

        #open files 
        fileObj = Path(fileName)
        if fileObj.is_file() == True:
            ds = xr.open_dataset(fileName)

        else: 
            print('File not found', fileName)
            continue



        #check if data values are realistic
        #lw
        ds['rlds_unrealistic_flag'] = (['time'], np.full((len(ds.rlds)), False))
        for i in range(0,len(ds.rlds)):
            if ds.rlds.values[i] > 700 or ds.rlds.values[i] < 10:
                ds.rlds_unrealistic_flag[i] = True 
                print('unrealistic rlds-values detected', day)

        #sw
        ds['rsds_unrealistic_flag'] = (['time'], np.full((len(ds.rsds)), False))
        for i in range(0, len(ds.rsds)):
            if ds.rsds.values[i] > 1500 or ds.rsds.values[i] < -5:
                ds.rsds_unrealistic_flag[i] = True 
                print('unrealistic rsds-values detected', day, print(ds.rsds.values[i]))



        #check for nan values
        ds['rlds_nan_flag'] = np.isnan(ds.rlds)
        if sum(ds['rlds_nan_flag'].data) > 0:
            print('nanrlds-Value detected', day)
            a = np.where(ds['rlds_nan_flag'] == True)
            print(len(a[0]), 'nan-rlds-Value detected', day, ds.time.values[a[0]])  
            #delete the timestamp with nan
            ds = ds.where(ds['rlds_nan_flag'] == False).dropna(how='all', dim='time')

        ds['rsds_nan_flag'] = np.isnan(ds.rsds)
        if sum(ds['rsds_nan_flag'].data) > 0:
            a = np.where(ds['rsds_nan_flag'] == True)
            print(len(a[0]), 'nan-rsds-Value detected', day, ds.time.values[a[0]])  
            #delete the timestamp with nan
            ds = ds.where(ds['rsds_nan_flag'] == False).dropna(how='all', dim='time')



        #check for gaps
        #create two lists with the timestamps, delete one object at the beginning and at the other list one at the end
        # substract the two list --> gaps length between one timestamp earlier and one later
        gap_end = ds.time.values
        gap_end = np.delete(gap_end, 0)
        gap_start = ds.time.values
        gap_start = np.delete(gap_start, len(gap_end)-1)

        gap_length = gap_end - gap_start
        
        #look where the gap is bigger than 5 seconds
        real_gaps = []
        #gives back the place of the gaps as a tuple
        real_gaps = np.where(gap_length > np.timedelta64(5000000000,'ns'))
        #print(real_gaps, len(real_gaps[0]), 'real_gaps')

        #if there are gaps, create a panda dataframe and save it
        #real_gaps is a tupel so len of the first vector
        if len(real_gaps[0]) > 0:

            gap_length_real = []
            gap_start_real = []
            gap_end_real = []
            missing_timestamps = []
            
            i = 0
            for i in range(0, len(real_gaps[0])):
                gap_length_real.append(gap_length[real_gaps[0][i]])
                gap_start_real.append(gap_start[real_gaps[0][i]])
                gap_end_real.append(gap_end[real_gaps[0][i]])
                missing_timestamps.append(gap_length[real_gaps[0][i]].astype('int64')/(5*10**9) -1)
                
            gap_start_real = gap_start_real +  np.timedelta64(5, 's')
            gap_end_real = gap_end_real -  np.timedelta64(5, 's')
            gap_length_real = gap_length_real - np.timedelta64(5, 's')
            
            
            data = {'start_of_gap': gap_start_real, 'end_of_gap': gap_end_real, 'length_of_gap': gap_length_real, 'missing_timestamps': missing_timestamps}
            df = pd.DataFrame(data=data)


            #checking if path/folder already exists
            MYDIR = (outpath + '/' + y)
            CHECK_FOLDER = os.path.isdir(MYDIR)

            #If folder doesn't exist, then create it.
            if not CHECK_FOLDER:
                os.makedirs(MYDIR)
                print("created folder : ", MYDIR)

            #create dataset and save it as an nc-File
            ds = df.to_xarray()
            
            path_and_name = MYDIR + '/' + outname + '_' +  str(y) + str(m) + str(d) + '.nc'
            ds.to_netcdf(path_and_name)  

            
            


# In[ ]:


#help(data_quality)
#help(merge_datasets)
startdate = '2018-01-01'
enddate = '2019-12-31'
merge_datasets(startdate, enddate, 'full_datasets', 'all_data')
data_quality(startdate, enddate, 'full_datasets', 'all_data', 'data_quality', 'gaps')

