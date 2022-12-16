#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import datetime as dt
from datetime import datetime, timedelta
from suntime import Sun, SunTimeException
from datetime import datetime, timezone
import numpy as np


def time_in_range(start, end, x):
    """Return true if x is in the range [start, end]"""
    if start <= end:
        return start <= x <= end
    else:
        return start <= x or x <= end   

def timesteps_light_juelich(startdate, enddate, path, name):
    '''
    What does it?
    -> In a timeperiod (between startdate and enddate) creates a list with 5 minutes timesteps, 
    during daylight (one hour after sunrise and one hour after sunrise),
    Saves it as a dataset, where Year, Month and Day are datavariables
    
    input:
        - startdate: <string>: e.g.: '2020-01-01' 
        - enddate:   <string>: e.g.: '2020-03-02'
        - path:      <string>: where to save the nc-file eg.: 'savings/'
        - name:      <string>: name of the nc-file: 'daylight.nc'
        
    ''' 
    #coordinates of Juelich
    latitude = 50.9224
    longitude = 6.3639
    sun = Sun(latitude, longitude)

    #period of time you want to have check for daylight
    all_dates = pd.date_range(start=startdate,end=enddate)

    light_flag = []
    date_string = []
    hour_string = []
    minute_string = []
    t_all = []

    for day in all_dates:

        #getting the time of sunrise +1 hour and the time of sunset - 1 hour. (result has information of timezone)
        sur_tzd = sun.get_local_sunrise_time(day) + timedelta(hours=1)
        sus_tzd = sun.get_local_sunset_time(day) - timedelta(hours=1)

        #delete timezonedesignator (otherwise it can't sompare if a certain time is within the range later (func time_in_range))
        sur = sur_tzd.replace(tzinfo=None)
        sus = sus_tzd.replace(tzinfo=None)

        #create an array with all timestamps (every five minutes) for the day
        next_day = day + timedelta(days=1)
        t = np.arange(day,next_day, timedelta(minutes=5)).astype(datetime)

        # check if the timestamps are between the sunrise and sunset
        for moment in t: 
            # check if the timestamps are between the sunrise and sunset
            start = sur
            end = sus
            light_flag.append(time_in_range(start, end, moment))

            #convert timstamps in strings (yyyymmdd, hh, mm)
            date_string.append(moment.strftime("%Y") + moment.strftime("%m") + moment.strftime("%d"))
            hour_string.append(moment.strftime('%H'))
            minute_string.append(moment.strftime('%M'))

        #create a list with all days
        t_all = np.concatenate((t_all, t))

    #create a dataframe
    d = {'date': date_string, 'hour': hour_string, 'minute': minute_string,  'light_flag': light_flag}
    df = pd.DataFrame(data=d)

    df.drop(df[df['light_flag'] == False].index, inplace = True)
    df = df.drop(columns=['light_flag'])

    #create dataset and save it as an nc-File
    ds = df.to_xarray()
    path_and_name = path + name
    ds.to_netcdf(path_and_name)  


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




