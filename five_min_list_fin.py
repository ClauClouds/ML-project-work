#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
import bisect
import matplotlib.pyplot as plt


# In[5]:




###################### all functions ###########################

def sunlight_juelich(day, ds, time2): # needed in merge_lw_sw_dataset
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



def merge_lw_sw_datasets(startdate, enddate, outpath, outname):
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



def data_quality(startdate, enddate, inpath, inname, outpath, outname, unrealistic_value):
    
    '''
    - you can choose a range of time, where it selects the datsets of the days where BOTH (LW and SW) Variables exists.
        (it uses the datasets which were created by the funct: merge_datasets)
    - checks the datasets for nan values (deletes the timesteps), unrealistic values (lw 5% & 95% percentil of year 2018 + 2019: 265 < 400) (0<sw), and for gaps (non existent timesteps)
    - creates a dataframe from all days with the following information: starttime, endtime, length, how many timesteps
    - creates a dataset with all days where no dataset was found for
    
     input:
        - startdate: <string>: e.g.: '2020-01-01' 
        - enddate:   <string>: e.g.: '2020-03-02'
        - inname:    <string>: e.g.: 'all_data'  automaticlly appends the date of the day and an .nc at the end
        - inpath:    <string>: where to save the nc-file eg.: 'full_datasets'
        - outpath:   <string>: where to save the csv-file eg.: 'data_quality'
        - outname:   <string>: name of the pkl file (date is added): eg.: 'gaps' --> gaps_20200101_20200302.pkl
        - unrealistic_value:<bool>  : True-> lw 5% & 95% percentil of year 2018 + 2019: 265 < 400,  0<sw
        
     output: 
        - saves gaps in a dataframe (starttime (the first missing), endtime (the last missing), length (counting one missing step as 5 sec), how many timesteps) as csv files (for each day)
    '''
    #Range of the datas
    all_dates = pd.date_range(start=startdate,end=enddate)
    ds_merge = None
    df_merge = None
    
    #coordinates of Juelich
    latitude = 50.908546
    longitude = 6.413536
    sun = Sun(latitude, longitude)
    
    # vektor for missing datasets 
    missing_dataset = np.full((len(all_dates)), 0)
    
    #loop backwards, because the df.merge, appends the dataframes on the top
    count = len(all_dates) - 1 
    for day in all_dates[::-1]:
        
        
        #get the day to string, to open the right file
        date = day.strftime('%Y%m%d') + '000000.nc'
        d = day.strftime('%d')
        m = day.strftime('%m')
        y = day.strftime('%Y')


        #Open Files for LW and shortwave data and put them in one dataset
        fileName = inpath +'/'+ y + '/'+ inname + date 
        
        #open files 
        fileObj = Path(fileName)
        if fileObj.is_file() == True:
            ds = xr.open_dataset(fileName)

        else: 
            print('File not found', fileName)
            #where datasets are missing put the value to -1
            missing_dataset[count] = -1
            count = count - 1
            continue
        

        if unrealistic_value == True:
            #check if data values are realistic
            #lw
            ds['rlds_unrealistic_flag'] = (['time'], np.full((len(ds.rlds)), False))
            for i in range(0,len(ds.rlds)):
                if ds.rlds.values[i] > 400 or ds.rlds.values[i] < 265:
                    ds.rlds_unrealistic_flag[i] = True 
                    #print('unrealistic rlds-values detected', day, ds.time.values[i], print(ds.rlds.values[i]))
            ds = ds.where(ds['rlds_unrealistic_flag'] == False).dropna(how='all', dim='time')


            #sw
            ds['rsds_unrealistic_flag'] = (['time'], np.full((len(ds.rsds)), False))
            for i in range(0, len(ds.rsds) - 1):

                if ds.rsds.values[i] > 3000 or ds.rsds.values[i] < 0:
                    ds.rsds_unrealistic_flag[i] = True 
                    #print('unrealistic rsds-values detected', day, ds.time.values[i], print(ds.rsds.values[i]))
            ds = ds.where(ds['rsds_unrealistic_flag'] == False).dropna(how='all', dim='time')
            



        #check for nan values
        ds['rlds_nan_flag'] = np.isnan(ds.rlds)
        if sum(ds['rlds_nan_flag'].data) > 0:
            #print('nanrlds-Value detected', day)
            a = np.where(ds['rlds_nan_flag'] == True)
            #print(len(a[0]), 'nan-rlds-Value detected', day, ds.time.values[a[0]])  
            #delete the timestamp with nan
            ds = ds.where(ds['rlds_nan_flag'] == False).dropna(how='all', dim='time')

        ds['rsds_nan_flag'] = np.isnan(ds.rsds)
        if sum(ds['rsds_nan_flag'].data) > 0:
            a = np.where(ds['rsds_nan_flag'] == True)
            #print(len(a[0]), 'nan-rsds-Value detected', day, ds.time.values[a[0]])  
            #delete the timestamp with nan
            ds = ds.where(ds['rsds_nan_flag'] == False).dropna(how='all', dim='time')



        #check for gaps
        #create two lists with the timestamps, delete one object at the beginning and at the other list one at the end
        # substract the two list --> gaps length between one timestamp earlier and one later
        if len(ds.time.values) == 0: #dataset is empty, because all entries have been deleted -> whole day is missing
            missing_dataset[count] = -1
            count = count - 1
            #print('file now is empty', day)
            continue
        count = count - 1
        
        gap_end = ds.time.values
        gap_end = np.delete(gap_end, 0)
        gap_start = ds.time.values
        gap_start = np.delete(gap_start, len(gap_end)-1)

        gap_length = gap_end - gap_start
        
        #look where the gap is bigger than 5 seconds
        real_gaps = []
        #gives back the place of the gaps as a tuple
        real_gaps = np.where(gap_length > np.timedelta64(5000000000,'ns'))

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
            
            #create dataframe and merge dataframes from all dates
            data = {'start_of_gap': gap_start_real, 'end_of_gap': gap_end_real, 'length_of_gap': gap_length_real, 'missing_timestamps': missing_timestamps}
            df = pd.DataFrame(data=data)
            
            #checking if path/folder already exists
            MYDIR = (outpath)
            CHECK_FOLDER = os.path.isdir(MYDIR)

            #If folder doesn't exist, then create it.
            if not CHECK_FOLDER:
                os.makedirs(MYDIR)
                print("created folder : ", MYDIR)

            #create dataset and save it as an nc-File
            ds_gaps = df.to_xarray()

            path_and_name = MYDIR + '/' + day.strftime('%Y') +'/' + outname + '_' + day.strftime('%Y%m%d') 
            ds_gaps.to_netcdf(path_and_name + '.nc') 
            df.to_pickle(path_and_name + '.pkl')

   
            if df_merge is None:
                df_merge = df
            else:
                df_merge = df.merge(df_merge, how= 'outer') 
            
    #checking if path/folder already exists
    MYDIR = (outpath)
    CHECK_FOLDER = os.path.isdir(MYDIR)

    #If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)
        
    path_and_name = MYDIR + '/' + outname  +  all_dates[0].strftime('%Y%m%d') + '_' + all_dates[len(all_dates)-1].strftime('%Y%m%d') 
        
    #save dataframe
    df_merge.to_pickle(path_and_name + '.pkl')
    df_merge.to_csv(path_and_name + '.csv')
    
    #save as dataset: create dataset and save it as an nc-File
    ds_merge = df_merge.to_xarray()
    ds_merge.to_netcdf(path_and_name + '.nc') 
    
    #save the dates where no dataset was available:
    data_missing = {'days' :all_dates, 'missing_dataset': missing_dataset}
    df_missing = pd.DataFrame(data=data_missing)
    df_missing.to_pickle(MYDIR + '/missing_datasets_'+  all_dates[0].strftime('%Y%m%d') + '_' + all_dates[len(all_dates)-1].strftime('%Y%m%d') + '.pkl')
    df_missing.to_csv(MYDIR + '/missing_datasets_'+  all_dates[0].strftime('%Y%m%d') + '_' + all_dates[len(all_dates)-1].strftime('%Y%m%d') + '.csv')

    

def gaps_per_day_in_a_year(year):
    
    '''
    What does it do?
        Plots the daily number of gaps in a year, when the plot is -1 then the dataset for the whole day is missing.
        Pre-condition: run the funct data_quality from 01-01 to 12-31 of the year 
    input: 
        - year: <string>  year 
    output: 
        - plot of the number of gaps'''
    
    startdate = year + '-01-01'
    enddate = year + '-12-31' 
    days = pd.date_range(start=startdate,end=enddate)


    # dataframes with gaps and missing days
    df_gaps = pd.read_pickle('data_quality/' + 'gaps_' + year + '0101_' + year + '1231.pkl')  
    df_miss_dataset = pd.read_pickle('data_quality/' + 'missing_datasets_'+year+'0101_'+year+'1231.pkl')  


    a = []
    df_gaps['date_of_gap'] = np.zeros(len(df_gaps['start_of_gap']))
    for i in range(0, len(df_gaps['start_of_gap'])):
        a.append(df_gaps['start_of_gap'][i].replace(hour=0, minute=0, second=0))
    df_gaps['date_of_gap'] = a


    number_of_gaps = df_miss_dataset['missing_dataset']
    for d in range(0, len(days)):
        for g in range(0,len(df_gaps['date_of_gap'])):
            if df_gaps['date_of_gap'][g] == days[d]:
                #print(df_gaps['date_of_gap'][g], days[d])
                number_of_gaps[d] = number_of_gaps[d] + 1


    data = {'day': days, 'number_of_gaps' : number_of_gaps}
    df = pd.DataFrame(data=data)   
    ax = df.plot(kind = 'scatter', x = 'day', y = 'number_of_gaps', figsize=(15,5))
    ax.set_title('Daily number of gaps in ' + year +' , -1 => whole day is missing')
    ax.set_ylabel('Number of gaps')
    ax.set_xlabel('Time')
    
    plt.show() 
    
    
        

####################################   five min steps functions   #######################################

#nearest and time_in_range are used in five_min_steps

def nearest(s, ts):
    # Given a presorted list of timestamps:  s = sorted(index)
    i = bisect.bisect_left(s, ts)
    return min(s[max(0, i-1): i+2], key=lambda t: abs(ts.astype('O') - t))

def time_in_range(start, end, x):
    """Return true if x is in the range [start, end]"""
    if start <= end:
        return start <= x <= end
    else:
        return start <= x or x <= end 

def five_min_steps(startdate, enddate, inname, outpath, outname):
    #coordinates of Juelich
    latitude = 50.9224
    longitude = 6.3639
    sun = Sun(latitude, longitude)

    #period of time you want tocheck for daylight
    all_dates = pd.date_range(start=startdate,end=enddate)

    light_flag = []
    gap_flag = []
    date_string = []
    hour_string = []
    minute_string = []
    t_all = []


    for day in all_dates:

        date = day.strftime('%Y%m%d')
        y= day.strftime('%Y')
        
        #check if the whole dataset is missing. If so, is goes to the next day
        index = []
        df_miss_dataset = pd.read_pickle('data_quality/' + 'missing_datasets_' + y + '0101_' + y + '1231.pkl')
        index = df_miss_dataset.index[df_miss_dataset['days'] == date].tolist()
        if df_miss_dataset['missing_dataset'][index[0]] == -1:
            continue
        
        path_and_name = outpath + '/' + y + '/' + outname + day.strftime('%Y%m%d') + '.nc'

        #getting the time of sunrise +1 hour and the time of sunset - 1 hour. (result has information of timezone)
        sur_tzd = sun.get_local_sunrise_time(day) + timedelta(hours=1)
        sus_tzd = sun.get_local_sunset_time(day) - timedelta(hours=1)
        #delete timezonedesignator (otherwise it can't sompare if a certain time is within the range later (func time_in_range))
        sur = sur_tzd.replace(tzinfo=None)
        sus = sus_tzd.replace(tzinfo=None)


        #create an array with all timestamps (every five minutes) for the day
        next_day = day + timedelta(days=1)
        t = np.arange(day,next_day, timedelta(minutes=5)).astype(datetime)
        date_string = []
        hour_string = []
        minute_string = []
        light_flag = []
        
        missing_flag = np.full(len(t), False)

        # check if the timestamps are between the sunrise and sunset
        start = sur
        end = sus
        for moment in t: 
            # check if the timestamps are between the sunrise and sunset
            light_flag.append(time_in_range(start, end, moment))

            #convert timstamps in strings (yyyymmdd, hh, mm)
            date_string.append(moment.strftime("%Y") + moment.strftime("%m") + moment.strftime("%d"))
            hour_string.append(moment.strftime('%H'))
            minute_string.append(moment.strftime('%M'))

    
        #open the dataset where the gaps are saved
        fileName = 'data_quality/' + y + '/' + inname + date +'.nc'
        fileObj = Path(str(fileName))
        if fileObj.is_file() == True:
            ds = xr.open_dataset(fileName)
        else: 
            print('File not found', fileName)
            
            #create a dataframe
            d = {'date': date_string, 'hour': hour_string, 'minute': minute_string,  'light_flag': light_flag, 'missing_flag': missing_flag}
            df = pd.DataFrame(data=d)

            df.drop(df[df['light_flag'] == False].index, inplace = True)
            df.drop(df[df['missing_flag'] == True].index, inplace = True)
            df = df.drop(columns=['light_flag'])
            df = df.drop(columns=['missing_flag'])

            #create dataset and save it as an nc-File
            ds = df.to_xarray()
            #checking if path/folder already exists
            MYDIR = (outpath + '/' + y + '/')
            CHECK_FOLDER = os.path.isdir(MYDIR)

            #If folder doesn't exist, then create it.
            if not CHECK_FOLDER:
                os.makedirs(MYDIR)
                print("created folder : ", MYDIR)

            #path_and_name = outpath + '/' + outname + day.strftime('%Y%m%d') + '.nc'
            ds.to_netcdf(path_and_name)  
            df.to_csv(outpath + '/' + y + '/' + outname + day.strftime('%Y%m%d') + '.csv')
            continue


        
        for i in range(0, len(ds.start_of_gap.values)):
            
            #search in t for nearest value before the start and end of the gap
            gap_start = ds.start_of_gap.values[i].astype('datetime64[s]')
            gap_end = ds.end_of_gap.values[i].astype('datetime64[s]')
            near_start = nearest(t, gap_start)
            near_end = nearest(t, gap_end)

            #get the position of the timestamps in t (to compare them t and the gap need the same timeformat)
            gap_start = pd.to_datetime(gap_start)
            gap_end = pd.to_datetime(gap_end)
            near_start_index = np.where(t == near_start)
            near_end_index = np.where(t == near_end)

            # if the closest timstamp in the 'five-min-list' (t), are the same for the start and end of the gap, 
            # the missing data, are only about one timestep in the five-min-list
            # than it can directly checked how many minutes are missing. If there are more than 2 minutes missing,
            # the timestamp will be deletet
            if near_start_index == near_end_index:
                
                if ds.length_of_gap.values[i] > 120*10**9:
                    missing_flag[near_end_index] = True
                    continue

            if gap_end - gap_start < timedelta(minutes = 2):
                continue


            #delete from behind to the front, so the indexes stay the same
            if gap_end - near_end > timedelta(minutes = 0):
                missing_flag[near_end_index] = True
            
            for index in range(near_end_index[0][0] - 1 , near_start_index[0][0], -1): #looping backwards 
                missing_flag[index] = True
                
            if gap_start - near_start < timedelta(minutes = - 0.5):
                missing_flag[near_start_index] = True

        #create a dataframe
        d = {'date': date_string, 'hour': hour_string, 'minute': minute_string,  'light_flag': light_flag, 'missing_flag': missing_flag}
        df = pd.DataFrame(data=d)

        df.drop(df[df['light_flag'] == False].index, inplace = True)
        df.drop(df[df['missing_flag'] == True].index, inplace = True)
        df = df.drop(columns=['light_flag'])
        df = df.drop(columns=['missing_flag'])

        #create dataset and save it as an nc-File
        ds = df.to_xarray()
        #checking if path/folder already exists
        MYDIR = (outpath + '/'+ y + '/' )
        CHECK_FOLDER = os.path.isdir(MYDIR)

        #If folder doesn't exist, then create it.
        if not CHECK_FOLDER:
            os.makedirs(MYDIR)
            print("created folder : ", MYDIR)

        #path_and_name = outpath + '/' + outname + day.strftime('%Y%m%d') + '.nc'
        ds.to_netcdf(path_and_name)  
        df.to_csv(outpath + '/' + y+ '/' + outname + day.strftime('%Y%m%d') + '.csv')
        
        
def merge_fivemin_datasets(startdate, enddate, inname, outpath, outname):
    '''
    merges the datasets of lw and sw in the given timerange to one big dataset.
    
    - outname: <string>: eg. 
    
    '''   
    inpath = 'savings/'
    
    
    ds_merge = None
    #Range of the datas
    all_dates = pd.date_range(start=startdate,end=enddate)

    for day in all_dates:
        y = day.strftime('%Y')
        
        
        date = day.strftime('%Y%m%d') + '.nc'
    
        
        #open dataset
        fileName =  inpath + y + '/' + inname + date 

        #open files 
        fileObj = Path(fileName)
        if fileObj.is_file() == True:
            ds = xr.open_dataset(fileName)
        else: 
            print('File not found', fileName)
            continue

        #print(ds_merge)
        if ds_merge is None:
            ds_merge = ds
        else:
            #ds_merge = ds.merge(ds_merge)
            ds_merge = xr.concat([ds_merge, ds], dim = 'index' )
        #except:
        #    print(f'doing ', day, 'could not merge! alarm!')
        
            
    #checking if path/folder already exists
    MYDIR = (outpath)
    CHECK_FOLDER = os.path.isdir(MYDIR)

    #If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)

    path_and_name = MYDIR + '/' + outname + '_' +  all_dates[0].strftime('%Y%m%d') + '_' + all_dates[len(all_dates)-1].strftime('%Y%m%d')  + '.nc'
    ds_merge.to_netcdf(path_and_name)  

    return ds_merge


# All three functions below create the Lists for with five min steps.
# The data quality function deletes:
# 
#      -  nan-Values
#      -  the 5% percentile of the lw radiation 
#      -  all sw values < 0     
# the gaps that this creates are saved. Also a list with missing datasets is created
# 
# The funct five_min_steps creates a list for every day with five minute steps from the sunrise (+1h) till the sunset (-1h), where the datasets from both lw & sw are existing.
# Then it looks where the gaps are, and deletes the five min timestamps, if the gap is lager than 2 minutes in the range of one timestamp. It assumes that the timestamp is in the middle. E.g. gap is from 12:31- 12:38 -> deletes 12:30 and 12:35
# 
# The func merge_fivemin_dataset merges all the daily datasets to a yearly dataset. 
# 
#      

# In[ ]:


#run data_quality function to work with the result from it

startdate = '2019-01-01'
enddate = '2019-12-31'
y = '2019'


inpath = 'full_datasets'
inname =  'all_data' 
outpath = 'data_quality'
outname = 'gaps_with_percentile'
unrealistic_value_percentile = True
#outname = 'gaps_'
#unrealistic_value_percentile = False

print('data_quality')
data_quality(startdate, enddate, inpath, inname, outpath, outname, unrealistic_value_percentile)


inname = 'gaps_with_percentile_'
outname = 'daylight_percentile' #appends an date.nc
#inname = 'gaps_'
#outname = 'daylight'
outpath = 'savings'
print('five_min')
five_min_steps(startdate, enddate,inname, outpath, outname)

inname = 'daylight_percentile'
outname = 'afive_min_steps_percentile'
#inname = 'daylight'
#outname = 'afive_min_steps'
print('merge_fivemin')
merge_fivemin_datasets(startdate, enddate, inname, 'savings', outname)


# In[ ]:




