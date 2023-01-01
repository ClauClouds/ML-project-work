#!/usr/bin/env python
# coding: utf-8

# In[4]:


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


def select_daylight_data_juelich(startdate, enddate, path):
    '''
    What does it?
    Step one from the instructions:
    
    Creates a script that reads the files sups_joy_pyr00 (SW data) and sups_joy_pyrg00 (LW data)
    from the folder data/hatpro/jue/hdcp2/radiation_hdcp2 and selects values at the requested satellite time stamps, 
    and saves them in a ncdf file
    
    -> In a timeperiod (between startdate and enddate) opens the daily datasets and selects the timstamps with sunlight.
    Saves it as datasets in a given path. It creates new folders for every year if they do not already exist.
    e.g.: 'path/2018/sups_joy_pyrg00_l1_rlds_v01_daylight20180102000000.nc'
    
    input:
        - startdate: <string>: e.g.: '2020-01-01' 
        - enddate:   <string>: e.g.: '2020-03-02'
        - path:      <string>: where to save the nc-file eg.: 'savings/'

        
    ''' 

    #Range of the datas
    all_dates = pd.date_range(start=startdate,end=enddate)

    #coordinates of Juelich
    latitude = 50.908546
    longitude = 6.413536
    sun = Sun(latitude, longitude)

    for day in all_dates:
        
        
        #getting the time of sunrise +1 hour and the time of sunset - 1 hour. (result has information of timezone)
        sur_tzd = sun.get_local_sunrise_time(day) + timedelta(hours=1)
        sus_tzd = sun.get_local_sunset_time(day) - timedelta(hours=1)

        #delete timezonedesignator (Data from sup_joys are in UTC)
        sur = sur_tzd.replace(tzinfo=None)
        sus = sus_tzd.replace(tzinfo=None)
        #print(sur, 'sunrise')
        

        
        #get the day to string, to open the right file
        y = day.strftime('%Y')
        m = day.strftime('%m')
        d = day.strftime('%d')
        date = y + m + d + '000000.nc'

        #Open Files for LW and shortwave data and put them in one dataset
        
        #starting with LW
        #getting the names from the data. Some names have ..._v00_... and some ..._v01...
        fileName_LW = glob.glob('/data/hatpro/jue/hdcp2/radiation_hdcp2/'+ y +'/sups_joy_pyrg00_l1_rlds_v*_'+ date )
        #if there is no file for the day it goes to the next day
        if not fileName_LW:
            print(date + ' no LW-file for this day found')
            continue
            
        #open files 
        fileObj_LW = Path(str(fileName_LW[0]))
        if fileObj_LW.is_file() == True:
            ds_LW = xr.open_dataset(str(fileName_LW[0]))
            
            #in a really few cases it didn't work out to select times (I guess something is wrong with the datasets). Nevertheless that the program doesent stop there is a try:...
            try:
                ds_LW_light = ds_LW.sel(time=slice(sur, sus))
            
            except:
                print(str(fileName_LW[0]),' something went wrong with selecting sunlight time')
           
            #checking if path/folder already exists
            MYDIR = (path + y)
            CHECK_FOLDER = os.path.isdir(MYDIR)

            # If folder doesn't exist, then create it.
            if not CHECK_FOLDER:
                os.makedirs(MYDIR)
                print("created folder : ", MYDIR)
                
            #saving it as a new dataset 
            save_LW = path + y + '/sups_joy_pyrg00_l1_rlds_v_daylight' + date
            ds_LW_light.to_netcdf(save_LW)  

        else: 
            print('File not found')
            print(fileName[0]_LW)

        #SW
        #getting the names from the data. Some names have ..._v00_... and some ..._v01...
        fileName_SW = glob.glob('/data/hatpro/jue/hdcp2/radiation_hdcp2/'+ y +'/sups_joy_pyr00_l1_rsds_v*_'+ date )
        #if there is no file for the day it goes to the next day
        if not fileName_SW:
            print(date + ' no SW-file for this day found')
            continue
            
        #open files 
        fileObj_SW = Path(str(fileName_SW[0]))
        if fileObj_SW.is_file() == True:
            ds_SW = xr.open_dataset(str(fileName_SW[0]))
            
            #in a really few cases it didn't work out to select times (I guess something is wrong with the datasets). Nevertheless that the program doesent stop there is a try:...
            try:
                ds_SW_light = ds_SW.sel(time=slice(sur, sus))
            
            except:
                print(str(fileName_LW[0]),' something went wrong with selecting sunlight time')
           
            #checking if path/folder already exists
            MYDIR = (path + y)
            CHECK_FOLDER = os.path.isdir(MYDIR)

            # If folder doesn't exist, then create it.
            if not CHECK_FOLDER:
                os.makedirs(MYDIR)
                print("created folder : ", MYDIR)
                
            #saving it as a new dataset 
            save_SW = path + y + '/sups_joy_pyr00_l1_rsds_v_daylight' + date
            ds_SW_light.to_netcdf(save_SW)  

        else: 
            print('File not found')
            print(fileName[0]) 

def select_daylight_solarpanel_data(startdate, enddate, path):
    '''
    What does it?
    Step two from the instructions
    
    -> In a timeperiod (between startdate and enddate) unzips the daily datasets, opens them and selects the timstamps with sunlight.
    Saves it as datasets in a given path. It creates new folders for every year if they do not already exist.
    e.g.: 'path/2016/details20161005.nc'
    
    input:
        - startdate: <string>: e.g.: '2020-01-01'  Only datas in Year 2016, 2017, 2020, 2021, 2022
        - enddate:   <string>: e.g.: '2020-03-02'
        - path:      <string>: where to save the nc-file eg.: 'savings1/'

        
    ''' 

    #Range of the datas
    all_dates = pd.date_range(start=startdate,end=enddate)

    #coordinates of Juelich
    latitude = 50.908546
    longitude = 6.413536
    sun = Sun(latitude, longitude)

    for day in all_dates:
        ds = xr.Dataset()
        ds_light = xr.Dataset()

        
        #getting the time of sunrise +1 hour and the time of sunset - 1 hour. (result has information of timezone)
        sur_tzd = sun.get_local_sunrise_time(day) + timedelta(hours=1)
        sus_tzd = sun.get_local_sunset_time(day) - timedelta(hours=1)

        #delete timezonedesignator (Data from sup_joys are in UTC)
        sur = sur_tzd.replace(tzinfo=None)
        sus = sus_tzd.replace(tzinfo=None)
        #print(sur, 'sunrise')
        

        
        #get the day to string, to open the right file
        y = day.strftime('%Y')
        m = day.strftime('%m')
        d = day.strftime('%d')
        date = y + m + d 


        #Unzipp the data and saving them in between
        fileName = '/data/obs/site/jue/pvm/l1/' + y + '/' + m + '/details' + date + '.nc.gz'
        fileObj = Path(fileName)
        if fileObj.is_file() == True:
            with gzip.open('/data/obs/site/jue/pvm/l1/' + y + '/' + m + '/details' + date + '.nc.gz', 'rb') as f_in:
                with open('date.nc', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
        else: 
            print('File not found: ', fileName)
            continue
       
    
        #open data
        if fileObj.is_file() == True:
            ds = xr.open_dataset('date.nc')
 
        else: 
            print(date, ' File not found: probably something went wrong with unzipping them')
            continue
            
        ds_light = ds.sel(time10s=slice(sur, sus), time60s=slice(sur, sus))
        

        
        #checking if path/folder already exists
        MYDIR = (path + y)
        CHECK_FOLDER = os.path.isdir(MYDIR)

        # If folder doesn't exist, then create it.
        if not CHECK_FOLDER:
            os.makedirs(MYDIR)
            print("created folder : ", MYDIR)

        #saving it as a new dataset 
        save = path + y + '/details' + date + '.nc'
        try:
            ds_light.to_netcdf(save)  
            #print(date, 'saved')
        except:
            print(date, 'saving didnt work')
        

    # Deleting the in between unzipped file
    file_path = 'date.nc'
    os.remove(file_path)



# In[5]:


select_daylight_data_juelich('2011-01-01', '2022-06-01', 'savings/')
  

select_daylight_solarpanel_data('2016-10-01', '2017-10-30', 'savings1/')


# In[6]:


xr.open_dataset('savings1/2016/details20161006' + '.nc')


# In[4]:


import pandas as pd
import datetime as dt
from datetime import datetime, timedelta
from suntime import Sun, SunTimeException
from datetime import datetime, timezone
import numpy as np
import glob

import os
import pandas as pd
import xarray as xr
from pathlib import Path



def new_Dataset(startdate, enddate, path, endname):
    '''
    What does it?
    Step 4 from the Instructions:
    
    4. Saving data in array dataset instead of pandas data frame
        1. Iterate on a list of files from input
        2. Read also SW and LW downwelling, as well as direct/diffuse SW radiation
    
    
    -> In a timeperiod (between startdate and enddate) opens the daily datasets, puts all data (LW, LW_error, SW, SW_error) 
    into one dataset and selects the timstamps with sunlight.
    Saves it as a dataset in a given name and path. It creates new folders for every year if they do not already exist.
    
    
    input:
        - startdate: <string>: e.g.: '2020-01-01' 
        - enddate:   <string>: e.g.: '2020-03-02'
        - endname:   <string>: e.g.: 'all_data'  automaticlly appends the date of the day and an .nc at the end
        - path:      <string>: where to save the nc-file eg.: 'savings/'

        
    ''' 



    #Range of the datas
    all_dates = pd.date_range(start=startdate,end=enddate)



    #coordinates of Juelich
    latitude = 50.908546
    longitude = 6.413536
    sun = Sun(latitude, longitude)

    for day in all_dates:

        #get the day to string, to open the right file
        y = day.strftime('%Y')
        m = day.strftime('%m')
        d = day.strftime('%d')
        date = y + m + d + '000000.nc'

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


        #getting the time of sunrise +1 hour and the time of sunset - 1 hour. (result has information of timezone)
        sur_tzd = sun.get_local_sunrise_time(day) + timedelta(hours=1)
        sus_tzd = sun.get_local_sunset_time(day) - timedelta(hours=1)

        #delete timezonedesignator (Data from sup_joys are in UTC)
        sur = sur_tzd.replace(tzinfo=None)
        sus = sus_tzd.replace(tzinfo=None)

        ds_light = ds.sel(time=slice(sur, sus))


        #checking if path/folder already exists
        MYDIR = (path + '/' + y)
        CHECK_FOLDER = os.path.isdir(MYDIR)

        # If folder doesn't exist, then create it.
        if not CHECK_FOLDER:
            os.makedirs(MYDIR)
            print("created folder : ", MYDIR)

        #saving it as a new dataset 
        save = path +'/'+ y + '/'+ endname + date 
        ds_light.to_netcdf(save)  


new_Dataset('2014-10-01', '2015-10-25', 'full_datasets', 'all_data')


# In[4]:


#how to open the dataset and read it to use it for the program

import xarray as xr
from pathlib import Path

startdate = '2020-01-01'
enddate = '2020-02-01'
path = 'savings/'
name = 'daylight'

#call function
timesteps_light_juelich(startdate, enddate, path, name)


fileName = path + name
fileObj = Path(fileName)
if fileObj.is_file() == True:
    ds = xr.open_dataset(fileName)
else: 
    print('File not found') 

dates = ds.date.data
hours = ds.hour.data
minutes = ds.minute.data

print(dates)
print(hours)
print(minutes)


# In[ ]:




