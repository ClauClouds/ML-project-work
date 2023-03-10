{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "9eac78a1",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "import datetime as dt\n",
    "from datetime import datetime, timedelta\n",
    "from suntime import Sun, SunTimeException\n",
    "from datetime import datetime, timezone\n",
    "import numpy as np\n",
    "import glob\n",
    "import os\n",
    "import xarray as xr\n",
    "from pathlib import Path\n",
    "import gzip\n",
    "import shutil\n",
    "import matplotlib.pyplot as plt"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "6b75eab7",
   "metadata": {},
   "outputs": [],
   "source": [
    "####################### all functions ###########################\n",
    "\n",
    "def sunlight_juelich(day, ds, time2):\n",
    "    '''\n",
    "    For the day it selects the timesteps with daylight (one hour after sunrise and one hour after sunrise)\n",
    "    !! time dimension from dataset needs to be called time or time10s and time60s!!\n",
    "    \n",
    "    time2: use for datasets from the solarpanel, for the rest put time2 = None'''\n",
    "\n",
    "    #coordinates of Juelich\n",
    "    latitude = 50.9224\n",
    "    longitude = 6.3639\n",
    "    sun = Sun(latitude, longitude)\n",
    "    \n",
    "    #getting the time of sunrise +1 hour and the time of sunset - 1 hour. (result has information of timezone)\n",
    "    sur_tzd = sun.get_local_sunrise_time(day) + timedelta(hours=1)\n",
    "    sus_tzd = sun.get_local_sunset_time(day) - timedelta(hours=1)\n",
    "\n",
    "    #delete timezonedesignator (Data from sup_joys are in UTC)\n",
    "    sur = sur_tzd.replace(tzinfo=None)\n",
    "    sus = sus_tzd.replace(tzinfo=None)\n",
    "    if time2 == None:\n",
    "        ds_light = ds.sel(time=slice(sur, sus))\n",
    "    else:\n",
    "        ds_light = ds.sel(time10s=slice(sur, sus), time60s=slice(sur, sus))\n",
    "   \n",
    "    return ds_light\n",
    "\n",
    "\n",
    "\n",
    "def merge_lw_sw_datasets(startdate, enddate, outpath, outname):\n",
    "    '''\n",
    "    What does it?\n",
    "    Step 4 from the Instructions:\n",
    "    \n",
    "    4. Saving data in array dataset instead of pandas data frame\n",
    "        1. Iterate on a list of files from input\n",
    "        2. Read also SW and LW downwelling, (as well as direct/diffuse SW radiation)\n",
    "    \n",
    "    \n",
    "    -> In a timeperiod (between startdate and enddate) opens the daily datasets, puts all data (LW, LW_error, SW, SW_error) \n",
    "    into one dataset and selects the timstamps with sunlight.\n",
    "    Saves it as a dataset in a given name and path. It creates new folders for every year if they do not already exist.\n",
    "    \n",
    "    \n",
    "    input:\n",
    "        - startdate: <string>: e.g.: '2020-01-01' \n",
    "        - enddate:   <string>: e.g.: '2020-03-02'\n",
    "        - outpath:   <string>: where to save the nc-file eg.: 'full_datasets'\n",
    "        - outname:   <string>: e.g.: 'all_data'  automaticlly appends the date of the day and an .nc at the end\n",
    "        \n",
    "        \n",
    "    output:\n",
    "        - saves a dataset (rlds, rlds_error, rsds, rsds_error), with the timesteps of the sunlight    \n",
    "    ''' \n",
    "    #Range of the datas\n",
    "    all_dates = pd.date_range(start=startdate,end=enddate)\n",
    "    \n",
    "\n",
    "    for day in all_dates:\n",
    "        \n",
    "        #get the day to string, to open the right file\n",
    "        date = day.strftime('%Y%m%d') + '000000.nc'\n",
    "        y = day.strftime('%Y')\n",
    "\n",
    "        #Open Files for LW and shortwave data and put them in one dataset\n",
    "         #getting the names from the data. Some names have ..._v00_... and some ..._v01...\n",
    "        fileName_lw = glob.glob('/data/hatpro/jue/hdcp2/radiation_hdcp2/'+ y +'/sups_joy_pyrg00_l1_rlds_v*_'+ date )\n",
    "        #if there is no file for the day it goes to the next day\n",
    "        if not fileName_lw:\n",
    "            print(date + ' no LW-file for this day found')\n",
    "            continue\n",
    "\n",
    "        #open files \n",
    "        fileObj_lw = Path(str(fileName_lw[0]))\n",
    "        if fileObj_lw.is_file() == True:\n",
    "            ds_lw = xr.open_dataset(str(fileName_lw[0]),\n",
    "                                   drop_variables = ['lon', 'lat'])\n",
    "\n",
    "        else: \n",
    "            print('File not found', fileName_lw[0])\n",
    "            continue\n",
    "\n",
    "        #getting the names from the data. Some names have ..._v00_... and some ..._v01...\n",
    "        fileName_sw = glob.glob('/data/hatpro/jue/hdcp2/radiation_hdcp2/'+ y +'/sups_joy_pyr00_l1_rsds_v*_'+ date )\n",
    "        #if there is no file for the day it goes to the next day\n",
    "        if not fileName_sw:\n",
    "            print(date + ' no SW-file for this day found')\n",
    "            continue\n",
    "\n",
    "        #open files \n",
    "        fileObj_sw = Path(str(fileName_sw[0]))\n",
    "        if fileObj_sw.is_file() == True:\n",
    "            ds_sw = xr.open_dataset(str(fileName_sw[0]))\n",
    "        else: \n",
    "            print('File not found', fileName_sw[0])\n",
    "            continue\n",
    "\n",
    "        #put the sw data to the lw data\n",
    "        ds_lw['rsds'] = ds_sw['rsds']\n",
    "        ds_lw['rsds_error'] = ds_sw['rsds_error']\n",
    "        ds = ds_lw\n",
    "        \n",
    "        #data with daylight\n",
    "        ds_light = sunlight_juelich(day, ds, None)\n",
    "\n",
    "        #checking if path/folder already exists\n",
    "        MYDIR = (outpath + '/' + y)\n",
    "        CHECK_FOLDER = os.path.isdir(MYDIR)\n",
    "\n",
    "        # If folder doesn't exist, then create it.\n",
    "        if not CHECK_FOLDER:\n",
    "            os.makedirs(MYDIR)\n",
    "            print(\"created folder : \", MYDIR)\n",
    "\n",
    "        #saving it as a new dataset \n",
    "        save = outpath +'/'+ y + '/'+ outname + date \n",
    "        ds_light.to_netcdf(save)  \n",
    "\n",
    " \n",
    "\n",
    "        \n",
    "def merge_all_datasets(startdate, enddate, outpath, outname):\n",
    "    '''\n",
    "    merges the datasets of lw and sw in the given timerange.\n",
    "    \n",
    "    - outname: <string>: eg. \n",
    "    \n",
    "    '''   \n",
    "    inpath = 'full_datasets'\n",
    "    inname = 'all_data'\n",
    "    \n",
    "    ds_merge = None\n",
    "    #Range of the datas\n",
    "    all_dates = pd.date_range(start=startdate,end=enddate)\n",
    "\n",
    "    for day in all_dates:\n",
    "        date = day.strftime('%Y%m%d') + '000000.nc'\n",
    "        y = day.strftime('%Y')\n",
    "        \n",
    "        #open dataset\n",
    "        fileName = inpath +'/'+ y + '/'+ inname + date \n",
    "\n",
    "        #open files \n",
    "        fileObj = Path(fileName)\n",
    "        if fileObj.is_file() == True:\n",
    "            ds = xr.open_dataset(fileName)\n",
    "        else: \n",
    "            print('File not found', fileName)\n",
    "            continue\n",
    "\n",
    "        #print(ds_merge)\n",
    "        if ds_merge is None:\n",
    "            ds_merge = ds\n",
    "        else:\n",
    "            ds_merge = ds.merge(ds_merge)\n",
    "        #except:\n",
    "        #    print(f'doing ', day, 'could not merge! alarm!')\n",
    "            \n",
    "    #checking if path/folder already exists\n",
    "    MYDIR = (outpath)\n",
    "    CHECK_FOLDER = os.path.isdir(MYDIR)\n",
    "\n",
    "    #If folder doesn't exist, then create it.\n",
    "    if not CHECK_FOLDER:\n",
    "        os.makedirs(MYDIR)\n",
    "        print(\"created folder : \", MYDIR)\n",
    "\n",
    "    path_and_name = MYDIR + '/' + outname + '_' +  all_dates[0].strftime('%Y%m%d') + '_' + all_dates[len(all_dates)-1].strftime('%Y%m%d')  + '.nc'\n",
    "    ds_merge.to_netcdf(path_and_name)  \n",
    "\n",
    "    return ds\n",
    "\n",
    "    \n",
    "\n",
    "\n",
    "def data_quality(startdate, enddate, inpath, inname, outpath, outname):\n",
    "    '''\n",
    "    - you can choose a range of time, where it selects the datsets of the days where BOTH (LW and SW) Variables exists.\n",
    "        (it uses the datasets which were created by the funct: merge_datasets)\n",
    "    - checks the datasets for nan values (deletes the timesteps), unrealistic values (10<lw<700, -5<sw<1500 ???), and for gaps (non existent timesteps)\n",
    "    - creates a dataframe from all days with the following information: starttime, endtime, length, how many timesteps\n",
    "    - creates a dataset with all days where no dataset was found for\n",
    "    \n",
    "     input:\n",
    "        - startdate: <string>: e.g.: '2020-01-01' \n",
    "        - enddate:   <string>: e.g.: '2020-03-02'\n",
    "        - inname:    <string>: e.g.: 'all_data'  automaticlly appends the date of the day and an .nc at the end\n",
    "        - inpath:    <string>: where to save the nc-file eg.: 'full_datasets'\n",
    "        - outpath:   <string>: where to save the csv-file eg.: 'data_quality'\n",
    "        - outname:   <string>: name of the pkl file (date is added): eg.: 'gaps' --> gaps_20200101_20200302.pkl\n",
    "        \n",
    "     output: \n",
    "        - saves gaps in a dataframe (starttime (the first missing), endtime (the last missing), length (counting one missing step as 5 sec), how many timesteps) as csv files (for each day)\n",
    "    '''\n",
    "    #Range of the datas\n",
    "    all_dates = pd.date_range(start=startdate,end=enddate)\n",
    "    ds_merge = None\n",
    "    df_merge = None\n",
    "    \n",
    "    #coordinates of Juelich\n",
    "    latitude = 50.908546\n",
    "    longitude = 6.413536\n",
    "    sun = Sun(latitude, longitude)\n",
    "    \n",
    "    # vektor for missing datasets \n",
    "    missing_dataset = np.full((len(all_dates)), 0)\n",
    "    \n",
    "    #loop backwards, because the df.merge, appends the dataframes on the top\n",
    "    count = len(all_dates) - 1 \n",
    "    for day in all_dates[::-1]:\n",
    "        \n",
    "\n",
    "        #get the day to string, to open the right file\n",
    "        date = day.strftime('%Y%m%d') + '000000.nc'\n",
    "        d = day.strftime('%d')\n",
    "        m = day.strftime('%m')\n",
    "        y = day.strftime('%Y')\n",
    "\n",
    "\n",
    "        #Open Files for LW and shortwave data and put them in one dataset\n",
    "         #getting the names from the data. Some names have ..._v00_... and some ..._v01...\n",
    "        #fileName ='full_datasets/' + y + '/all_data' + date\n",
    "        fileName = inpath +'/'+ y + '/'+ inname + date \n",
    "        \n",
    "        #open files \n",
    "        fileObj = Path(fileName)\n",
    "        if fileObj.is_file() == True:\n",
    "            ds = xr.open_dataset(fileName)\n",
    "\n",
    "        else: \n",
    "            print('File not found', fileName)\n",
    "            #where datasets are missing put the value to -1\n",
    "            missing_dataset[count] = -1\n",
    "            count = count - 1\n",
    "            continue\n",
    "        count = count - 1\n",
    "\n",
    "\n",
    "\n",
    "        #check if data values are realistic\n",
    "        #lw\n",
    "        ds['rlds_unrealistic_flag'] = (['time'], np.full((len(ds.rlds)), False))\n",
    "        for i in range(0,len(ds.rlds)):\n",
    "            if ds.rlds.values[i] > 700 or ds.rlds.values[i] < 10:\n",
    "                ds.rlds_unrealistic_flag[i] = True \n",
    "                print('unrealistic rlds-values detected', day)\n",
    "\n",
    "        #sw\n",
    "        ds['rsds_unrealistic_flag'] = (['time'], np.full((len(ds.rsds)), False))\n",
    "        for i in range(0, len(ds.rsds)):\n",
    "            if ds.rsds.values[i] > 1500 or ds.rsds.values[i] < -5:\n",
    "                ds.rsds_unrealistic_flag[i] = True \n",
    "                print('unrealistic rsds-values detected', day, print(ds.rsds.values[i]))\n",
    "\n",
    "\n",
    "\n",
    "        #check for nan values\n",
    "        ds['rlds_nan_flag'] = np.isnan(ds.rlds)\n",
    "        if sum(ds['rlds_nan_flag'].data) > 0:\n",
    "            print('nanrlds-Value detected', day)\n",
    "            a = np.where(ds['rlds_nan_flag'] == True)\n",
    "            print(len(a[0]), 'nan-rlds-Value detected', day, ds.time.values[a[0]])  \n",
    "            #delete the timestamp with nan\n",
    "            ds = ds.where(ds['rlds_nan_flag'] == False).dropna(how='all', dim='time')\n",
    "\n",
    "        ds['rsds_nan_flag'] = np.isnan(ds.rsds)\n",
    "        if sum(ds['rsds_nan_flag'].data) > 0:\n",
    "            a = np.where(ds['rsds_nan_flag'] == True)\n",
    "            print(len(a[0]), 'nan-rsds-Value detected', day, ds.time.values[a[0]])  \n",
    "            #delete the timestamp with nan\n",
    "            ds = ds.where(ds['rsds_nan_flag'] == False).dropna(how='all', dim='time')\n",
    "\n",
    "\n",
    "\n",
    "        #check for gaps\n",
    "        #create two lists with the timestamps, delete one object at the beginning and at the other list one at the end\n",
    "        # substract the two list --> gaps length between one timestamp earlier and one later\n",
    "        gap_end = ds.time.values\n",
    "        gap_end = np.delete(gap_end, 0)\n",
    "        gap_start = ds.time.values\n",
    "        gap_start = np.delete(gap_start, len(gap_end)-1)\n",
    "\n",
    "        gap_length = gap_end - gap_start\n",
    "        \n",
    "        #look where the gap is bigger than 5 seconds\n",
    "        real_gaps = []\n",
    "        #gives back the place of the gaps as a tuple\n",
    "        real_gaps = np.where(gap_length > np.timedelta64(5000000000,'ns'))\n",
    "        #print(real_gaps, len(real_gaps[0]), 'real_gaps')\n",
    "\n",
    "        #if there are gaps, create a panda dataframe and save it\n",
    "        #real_gaps is a tupel so len of the first vector\n",
    "        if len(real_gaps[0]) > 0:\n",
    "\n",
    "            gap_length_real = []\n",
    "            gap_start_real = []\n",
    "            gap_end_real = []\n",
    "            missing_timestamps = []\n",
    "            \n",
    "            i = 0\n",
    "            for i in range(0, len(real_gaps[0])):\n",
    "                gap_length_real.append(gap_length[real_gaps[0][i]])\n",
    "                gap_start_real.append(gap_start[real_gaps[0][i]])\n",
    "                gap_end_real.append(gap_end[real_gaps[0][i]])\n",
    "                missing_timestamps.append(gap_length[real_gaps[0][i]].astype('int64')/(5*10**9) -1)\n",
    "                \n",
    "            gap_start_real = gap_start_real +  np.timedelta64(5, 's')\n",
    "            gap_end_real = gap_end_real -  np.timedelta64(5, 's')\n",
    "            gap_length_real = gap_length_real - np.timedelta64(5, 's')\n",
    "            \n",
    "            #create dataframe and merge dataframes from all dates\n",
    "            data = {'start_of_gap': gap_start_real, 'end_of_gap': gap_end_real, 'length_of_gap': gap_length_real, 'missing_timestamps': missing_timestamps}\n",
    "            df = pd.DataFrame(data=data)\n",
    "\n",
    "   \n",
    "            if df_merge is None:\n",
    "                df_merge = df\n",
    "            else:\n",
    "                df_merge = df.merge(df_merge, how= 'outer')\n",
    "            #print(df_merge)\n",
    "       \n",
    "            \n",
    "    #checking if path/folder already exists\n",
    "    MYDIR = (outpath)\n",
    "    CHECK_FOLDER = os.path.isdir(MYDIR)\n",
    "\n",
    "    #If folder doesn't exist, then create it.\n",
    "    if not CHECK_FOLDER:\n",
    "        os.makedirs(MYDIR)\n",
    "        print(\"created folder : \", MYDIR)\n",
    "        \n",
    "    path_and_name = MYDIR + '/' + outname + '_' +  all_dates[0].strftime('%Y%m%d') + '_' + all_dates[len(all_dates)-1].strftime('%Y%m%d') \n",
    "        \n",
    "    #save dataframe\n",
    "    df_merge.to_pickle(path_and_name + '.pkl')\n",
    "    \n",
    "    #save as dataset: create dataset and save it as an nc-File\n",
    "    ds_merge = df_merge.to_xarray()\n",
    "    ds_merge.to_netcdf(path_and_name + '.nc') \n",
    "    \n",
    "    #save the dates where no dataset was available:\n",
    "    data_missing = {'days' :all_dates, 'missing_dataset': missing_dataset}\n",
    "    df_missing = pd.DataFrame(data=data_missing)\n",
    "    df_missing.to_pickle(MYDIR + '/missing_datasets_'+  all_dates[0].strftime('%Y%m%d') + '_' + all_dates[len(all_dates)-1].strftime('%Y%m%d') + '.pkl')\n",
    "\n",
    "    \n",
    "def gaps_per_day_in_a_year(year):\n",
    "    \n",
    "    '''\n",
    "    What does it do?\n",
    "        Plots the daily number of gaps in a year, when the plot is -1 then the dataset for the whole day is missing.\n",
    "        Pre-condition: run the funct data_quality from 01-01 to 12-31 of the year \n",
    "    input: \n",
    "        - year: <string>  year \n",
    "    output: \n",
    "        - plot of the number of gaps'''\n",
    "    \n",
    "    startdate = year + '-01-01'\n",
    "    enddate = year + '-12-31' \n",
    "    days = pd.date_range(start=startdate,end=enddate)\n",
    "\n",
    "\n",
    "    # dataframes with gaps and missing days\n",
    "    df_gaps = pd.read_pickle('data_quality/' + 'gaps_' + year + '0101_' + year + '1231.pkl')  \n",
    "    df_miss_dataset = pd.read_pickle('data_quality/' + 'missing_datasets_'+year+'0101_'+year+'1231.pkl')  \n",
    "\n",
    "\n",
    "    a = []\n",
    "    df_gaps['date_of_gap'] = np.zeros(len(df_gaps['start_of_gap']))\n",
    "    for i in range(0, len(df_gaps['start_of_gap'])):\n",
    "        a.append(df_gaps['start_of_gap'][i].replace(hour=0, minute=0, second=0))\n",
    "    df_gaps['date_of_gap'] = a\n",
    "\n",
    "\n",
    "    number_of_gaps = df_miss_dataset['missing_dataset']\n",
    "    for d in range(0, len(days)):\n",
    "        for g in range(0,len(df_gaps['date_of_gap'])):\n",
    "            if df_gaps['date_of_gap'][g] == days[d]:\n",
    "                #print(df_gaps['date_of_gap'][g], days[d])\n",
    "                number_of_gaps[d] = number_of_gaps[d] + 1\n",
    "\n",
    "\n",
    "    data = {'day': days, 'number_of_gaps' : number_of_gaps}\n",
    "    df = pd.DataFrame(data=data)   \n",
    "    ax = df.plot(kind = 'scatter', x = 'day', y = 'number_of_gaps', figsize=(15,5))\n",
    "    ax.set_title('Daily number of gaps in ' + year +' , -1 => whole day is missing')\n",
    "    ax.set_ylabel('Number of gaps')\n",
    "    ax.set_xlabel('Time')\n",
    "    \n",
    "    plt.show() "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "a4a399e7",
   "metadata": {
    "collapsed": true
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "3 nan-rsds-Value detected 2019-07-25 00:00:00 ['2019-07-25T20:06:13.000000000' '2019-07-25T20:06:18.000000000'\n",
      " '2019-07-25T20:06:33.000000000']\n",
      "82 nan-rsds-Value detected 2019-07-22 00:00:00 ['2019-07-22T20:02:28.000000000' '2019-07-22T20:02:33.000000000'\n",
      " '2019-07-22T20:02:38.000000000' '2019-07-22T20:02:43.000000000'\n",
      " '2019-07-22T20:02:48.000000000' '2019-07-22T20:02:53.000000000'\n",
      " '2019-07-22T20:02:58.000000000' '2019-07-22T20:03:03.000000000'\n",
      " '2019-07-22T20:03:08.000000000' '2019-07-22T20:03:13.000000000'\n",
      " '2019-07-22T20:03:18.000000000' '2019-07-22T20:03:23.000000000'\n",
      " '2019-07-22T20:03:28.000000000' '2019-07-22T20:03:33.000000000'\n",
      " '2019-07-22T20:03:38.000000000' '2019-07-22T20:03:43.000000000'\n",
      " '2019-07-22T20:03:48.000000000' '2019-07-22T20:03:53.000000000'\n",
      " '2019-07-22T20:03:58.000000000' '2019-07-22T20:04:03.000000000'\n",
      " '2019-07-22T20:04:08.000000000' '2019-07-22T20:04:13.000000000'\n",
      " '2019-07-22T20:04:18.000000000' '2019-07-22T20:04:23.000000000'\n",
      " '2019-07-22T20:04:28.000000000' '2019-07-22T20:04:33.000000000'\n",
      " '2019-07-22T20:04:38.000000000' '2019-07-22T20:18:13.000000000'\n",
      " '2019-07-22T20:18:18.000000000' '2019-07-22T20:18:23.000000000'\n",
      " '2019-07-22T20:18:28.000000000' '2019-07-22T20:18:33.000000000'\n",
      " '2019-07-22T20:18:38.000000000' '2019-07-22T20:18:43.000000000'\n",
      " '2019-07-22T20:18:48.000000000' '2019-07-22T20:18:53.000000000'\n",
      " '2019-07-22T20:18:58.000000000' '2019-07-22T20:19:03.000000000'\n",
      " '2019-07-22T20:19:08.000000000' '2019-07-22T20:19:13.000000000'\n",
      " '2019-07-22T20:19:18.000000000' '2019-07-22T20:19:23.000000000'\n",
      " '2019-07-22T20:19:28.000000000' '2019-07-22T20:19:33.000000000'\n",
      " '2019-07-22T20:19:38.000000000' '2019-07-22T20:19:43.000000000'\n",
      " '2019-07-22T20:19:48.000000000' '2019-07-22T20:19:53.000000000'\n",
      " '2019-07-22T20:19:58.000000000' '2019-07-22T20:20:03.000000000'\n",
      " '2019-07-22T20:20:08.000000000' '2019-07-22T20:20:13.000000000'\n",
      " '2019-07-22T20:20:18.000000000' '2019-07-22T20:20:23.000000000'\n",
      " '2019-07-22T20:20:28.000000000' '2019-07-22T20:20:33.000000000'\n",
      " '2019-07-22T20:20:38.000000000' '2019-07-22T20:20:43.000000000'\n",
      " '2019-07-22T20:20:48.000000000' '2019-07-22T20:20:53.000000000'\n",
      " '2019-07-22T20:20:58.000000000' '2019-07-22T20:21:03.000000000'\n",
      " '2019-07-22T20:21:08.000000000' '2019-07-22T20:21:13.000000000'\n",
      " '2019-07-22T20:21:18.000000000' '2019-07-22T20:21:23.000000000'\n",
      " '2019-07-22T20:21:28.000000000' '2019-07-22T20:21:33.000000000'\n",
      " '2019-07-22T20:21:38.000000000' '2019-07-22T20:21:43.000000000'\n",
      " '2019-07-22T20:21:48.000000000' '2019-07-22T20:21:53.000000000'\n",
      " '2019-07-22T20:21:58.000000000' '2019-07-22T20:22:03.000000000'\n",
      " '2019-07-22T20:22:08.000000000' '2019-07-22T20:22:13.000000000'\n",
      " '2019-07-22T20:22:18.000000000' '2019-07-22T20:22:23.000000000'\n",
      " '2019-07-22T20:22:28.000000000' '2019-07-22T20:22:33.000000000'\n",
      " '2019-07-22T20:22:38.000000000' '2019-07-22T20:22:43.000000000']\n",
      "6 nan-rsds-Value detected 2019-06-28 00:00:00 ['2019-06-28T09:36:08.000000000' '2019-06-28T09:36:13.000000000'\n",
      " '2019-06-28T09:36:18.000000000' '2019-06-28T09:36:23.000000000'\n",
      " '2019-06-28T09:36:28.000000000' '2019-06-28T09:36:33.000000000']\n",
      "9 nan-rsds-Value detected 2019-06-27 00:00:00 ['2019-06-27T16:38:13.000000000' '2019-06-27T16:38:18.000000000'\n",
      " '2019-06-27T16:38:23.000000000' '2019-06-27T16:38:28.000000000'\n",
      " '2019-06-27T16:38:33.000000000' '2019-06-27T16:38:38.000000000'\n",
      " '2019-06-27T16:38:43.000000000' '2019-06-27T16:38:48.000000000'\n",
      " '2019-06-27T16:38:53.000000000']\n",
      "File not found full_datasets/2019/all_data20190506000000.nc\n"
     ]
    }
   ],
   "source": [
    "#run data_quality function to work with the result from it\n",
    "\n",
    "startdate = '2019-01-01'\n",
    "enddate = '2019-12-31'\n",
    "inpath = 'full_datasets'\n",
    "inname =  'all_data' \n",
    "outpath = 'data_quality'\n",
    "outname = 'gaps'\n",
    "data_quality(startdate, enddate, inpath, inname, outpath, outname)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "47195bac",
   "metadata": {},
   "source": [
    "### Plot the lw and sw in a historgram ###\n",
    "all data, without filters"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "id": "a51ed246",
   "metadata": {
    "collapsed": true
   },
   "outputs": [
    {
     "ename": "KeyboardInterrupt",
     "evalue": "",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mKeyboardInterrupt\u001b[0m                         Traceback (most recent call last)",
      "\u001b[0;32m/tmp/ipykernel_2887736/945240113.py\u001b[0m in \u001b[0;36m<module>\u001b[0;34m\u001b[0m\n\u001b[1;32m      2\u001b[0m \u001b[0menddate\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m'2019-12-31'\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m      3\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m----> 4\u001b[0;31m \u001b[0mmerge_all_datasets\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mstartdate\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0menddate\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'merged_dataset'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'merged'\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m",
      "\u001b[0;32m/tmp/ipykernel_2887736/2641684810.py\u001b[0m in \u001b[0;36mmerge_all_datasets\u001b[0;34m(startdate, enddate, outpath, outname)\u001b[0m\n\u001b[1;32m    154\u001b[0m             \u001b[0mds_merge\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mds\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    155\u001b[0m         \u001b[0;32melse\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 156\u001b[0;31m             \u001b[0mds_merge\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mmerge\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mds_merge\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    157\u001b[0m         \u001b[0;31m#except:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    158\u001b[0m         \u001b[0;31m#    print(f'doing ', day, 'could not merge! alarm!')\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/dataset.py\u001b[0m in \u001b[0;36mmerge\u001b[0;34m(self, other, overwrite_vars, compat, join, fill_value, combine_attrs)\u001b[0m\n\u001b[1;32m   5027\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   5028\u001b[0m         \u001b[0mother\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mother\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mto_dataset\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;32mif\u001b[0m \u001b[0misinstance\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mother\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mDataArray\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;32melse\u001b[0m \u001b[0mother\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m-> 5029\u001b[0;31m         merge_result = dataset_merge_method(\n\u001b[0m\u001b[1;32m   5030\u001b[0m             \u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   5031\u001b[0m             \u001b[0mother\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/merge.py\u001b[0m in \u001b[0;36mdataset_merge_method\u001b[0;34m(dataset, other, overwrite_vars, compat, join, fill_value, combine_attrs)\u001b[0m\n\u001b[1;32m   1069\u001b[0m         \u001b[0mpriority_arg\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m2\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   1070\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m-> 1071\u001b[0;31m     return merge_core(\n\u001b[0m\u001b[1;32m   1072\u001b[0m         \u001b[0mobjs\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   1073\u001b[0m         \u001b[0mcompat\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/merge.py\u001b[0m in \u001b[0;36mmerge_core\u001b[0;34m(objects, compat, join, combine_attrs, priority_arg, explicit_coords, indexes, fill_value)\u001b[0m\n\u001b[1;32m    755\u001b[0m     \u001b[0mcollected\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mcollect_variables_and_indexes\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0maligned\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mindexes\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mindexes\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    756\u001b[0m     \u001b[0mprioritized\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0m_get_priority_vars_and_indexes\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0maligned\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mpriority_arg\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mcompat\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mcompat\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 757\u001b[0;31m     variables, out_indexes = merge_collected(\n\u001b[0m\u001b[1;32m    758\u001b[0m         \u001b[0mcollected\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mprioritized\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mcompat\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mcompat\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mcombine_attrs\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mcombine_attrs\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    759\u001b[0m     )\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/merge.py\u001b[0m in \u001b[0;36mmerge_collected\u001b[0;34m(grouped, prioritized, compat, combine_attrs, equals)\u001b[0m\n\u001b[1;32m    300\u001b[0m                 \u001b[0mvariables\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m[\u001b[0m\u001b[0mvariable\u001b[0m \u001b[0;32mfor\u001b[0m \u001b[0mvariable\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0m_\u001b[0m \u001b[0;32min\u001b[0m \u001b[0melements_list\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    301\u001b[0m                 \u001b[0;32mtry\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 302\u001b[0;31m                     merged_vars[name] = unique_variable(\n\u001b[0m\u001b[1;32m    303\u001b[0m                         \u001b[0mname\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mvariables\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mcompat\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mequals\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mget\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mname\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;32mNone\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    304\u001b[0m                     )\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/merge.py\u001b[0m in \u001b[0;36munique_variable\u001b[0;34m(name, variables, compat, equals)\u001b[0m\n\u001b[1;32m    161\u001b[0m     \u001b[0;32mif\u001b[0m \u001b[0mcombine_method\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    162\u001b[0m         \u001b[0;32mfor\u001b[0m \u001b[0mvar\u001b[0m \u001b[0;32min\u001b[0m \u001b[0mvariables\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0;36m1\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 163\u001b[0;31m             \u001b[0mout\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mgetattr\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mout\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mcombine_method\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mvar\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    164\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    165\u001b[0m     \u001b[0;32mreturn\u001b[0m \u001b[0mout\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/variable.py\u001b[0m in \u001b[0;36mfillna\u001b[0;34m(self, value)\u001b[0m\n\u001b[1;32m   1867\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   1868\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0mfillna\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mvalue\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m-> 1869\u001b[0;31m         \u001b[0;32mreturn\u001b[0m \u001b[0mops\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mfillna\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mvalue\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m   1870\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   1871\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0mwhere\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mcond\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mother\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mdtypes\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mNA\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/ops.py\u001b[0m in \u001b[0;36mfillna\u001b[0;34m(data, other, join, dataset_join)\u001b[0m\n\u001b[1;32m    144\u001b[0m     \u001b[0;32mfrom\u001b[0m \u001b[0;34m.\u001b[0m\u001b[0mcomputation\u001b[0m \u001b[0;32mimport\u001b[0m \u001b[0mapply_ufunc\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    145\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 146\u001b[0;31m     return apply_ufunc(\n\u001b[0m\u001b[1;32m    147\u001b[0m         \u001b[0mduck_array_ops\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mfillna\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    148\u001b[0m         \u001b[0mdata\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/computation.py\u001b[0m in \u001b[0;36mapply_ufunc\u001b[0;34m(func, input_core_dims, output_core_dims, exclude_dims, vectorize, join, dataset_join, dataset_fill_value, keep_attrs, kwargs, dask, output_dtypes, output_sizes, meta, dask_gufunc_kwargs, *args)\u001b[0m\n\u001b[1;32m   1212\u001b[0m     \u001b[0;31m# feed Variables directly through apply_variable_ufunc\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   1213\u001b[0m     \u001b[0;32melif\u001b[0m \u001b[0many\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0misinstance\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0ma\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mVariable\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;32mfor\u001b[0m \u001b[0ma\u001b[0m \u001b[0;32min\u001b[0m \u001b[0margs\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m-> 1214\u001b[0;31m         \u001b[0;32mreturn\u001b[0m \u001b[0mvariables_vfunc\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m*\u001b[0m\u001b[0margs\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m   1215\u001b[0m     \u001b[0;32melse\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   1216\u001b[0m         \u001b[0;31m# feed anything else through apply_array_ufunc\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/computation.py\u001b[0m in \u001b[0;36mapply_variable_ufunc\u001b[0;34m(func, signature, exclude_dims, dask, output_dtypes, vectorize, keep_attrs, dask_gufunc_kwargs, *args)\u001b[0m\n\u001b[1;32m    769\u001b[0m             )\n\u001b[1;32m    770\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 771\u001b[0;31m     \u001b[0mresult_data\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mfunc\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m*\u001b[0m\u001b[0minput_data\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    772\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    773\u001b[0m     \u001b[0;32mif\u001b[0m \u001b[0msignature\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mnum_outputs\u001b[0m \u001b[0;34m==\u001b[0m \u001b[0;36m1\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/duck_array_ops.py\u001b[0m in \u001b[0;36mfillna\u001b[0;34m(data, other)\u001b[0m\n\u001b[1;32m    307\u001b[0m     \u001b[0;31m# correct unit\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    308\u001b[0m     \u001b[0;31m# TODO: revert after https://github.com/hgrecco/pint/issues/1019 is fixed\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 309\u001b[0;31m     \u001b[0;32mreturn\u001b[0m \u001b[0mwhere\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mnotnull\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mdata\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mdata\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mother\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    310\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    311\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/duck_array_ops.py\u001b[0m in \u001b[0;36mwhere\u001b[0;34m(condition, x, y)\u001b[0m\n\u001b[1;32m    294\u001b[0m     \u001b[0;34m\"\"\"Three argument where() with better dtype promotion rules.\"\"\"\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    295\u001b[0m     \u001b[0mxp\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mget_array_namespace\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mcondition\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 296\u001b[0;31m     \u001b[0;32mreturn\u001b[0m \u001b[0mxp\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mwhere\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mcondition\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m*\u001b[0m\u001b[0mas_shared_dtype\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mx\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0my\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mxp\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mxp\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    297\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    298\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/duck_array_ops.py\u001b[0m in \u001b[0;36mas_shared_dtype\u001b[0;34m(scalars_or_arrays, xp)\u001b[0m\n\u001b[1;32m    194\u001b[0m \u001b[0;32mdef\u001b[0m \u001b[0mas_shared_dtype\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mscalars_or_arrays\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mxp\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mnp\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    195\u001b[0m     \u001b[0;34m\"\"\"Cast a arrays to a shared dtype using xarray's type promotion rules.\"\"\"\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 196\u001b[0;31m     \u001b[0;32mif\u001b[0m \u001b[0many\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0misinstance\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mx\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0marray_type\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m\"cupy\"\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;32mfor\u001b[0m \u001b[0mx\u001b[0m \u001b[0;32min\u001b[0m \u001b[0mscalars_or_arrays\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    197\u001b[0m         \u001b[0;32mimport\u001b[0m \u001b[0mcupy\u001b[0m \u001b[0;32mas\u001b[0m \u001b[0mcp\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    198\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/duck_array_ops.py\u001b[0m in \u001b[0;36m<genexpr>\u001b[0;34m(.0)\u001b[0m\n\u001b[1;32m    194\u001b[0m \u001b[0;32mdef\u001b[0m \u001b[0mas_shared_dtype\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mscalars_or_arrays\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mxp\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mnp\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    195\u001b[0m     \u001b[0;34m\"\"\"Cast a arrays to a shared dtype using xarray's type promotion rules.\"\"\"\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 196\u001b[0;31m     \u001b[0;32mif\u001b[0m \u001b[0many\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0misinstance\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mx\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0marray_type\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m\"cupy\"\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;32mfor\u001b[0m \u001b[0mx\u001b[0m \u001b[0;32min\u001b[0m \u001b[0mscalars_or_arrays\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    197\u001b[0m         \u001b[0;32mimport\u001b[0m \u001b[0mcupy\u001b[0m \u001b[0;32mas\u001b[0m \u001b[0mcp\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    198\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/pycompat.py\u001b[0m in \u001b[0;36marray_type\u001b[0;34m(mod)\u001b[0m\n\u001b[1;32m     62\u001b[0m \u001b[0;32mdef\u001b[0m \u001b[0marray_type\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mmod\u001b[0m\u001b[0;34m:\u001b[0m \u001b[0mModType\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;34m->\u001b[0m \u001b[0mDuckArrayTypes\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     63\u001b[0m     \u001b[0;34m\"\"\"Quick wrapper to get the array class of the module.\"\"\"\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m---> 64\u001b[0;31m     \u001b[0;32mreturn\u001b[0m \u001b[0mDuckArrayModule\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mmod\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mtype\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m     65\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     66\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/pycompat.py\u001b[0m in \u001b[0;36m__init__\u001b[0;34m(self, mod)\u001b[0m\n\u001b[1;32m     35\u001b[0m         \u001b[0mduck_array_type\u001b[0m\u001b[0;34m:\u001b[0m \u001b[0mDuckArrayTypes\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     36\u001b[0m         \u001b[0;32mtry\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m---> 37\u001b[0;31m             \u001b[0mduck_array_module\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mimport_module\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mmod\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m     38\u001b[0m             \u001b[0mduck_array_version\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mVersion\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mduck_array_module\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m__version__\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     39\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/__init__.py\u001b[0m in \u001b[0;36mimport_module\u001b[0;34m(name, package)\u001b[0m\n\u001b[1;32m    124\u001b[0m                 \u001b[0;32mbreak\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    125\u001b[0m             \u001b[0mlevel\u001b[0m \u001b[0;34m+=\u001b[0m \u001b[0;36m1\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 126\u001b[0;31m     \u001b[0;32mreturn\u001b[0m \u001b[0m_bootstrap\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_gcd_import\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mname\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mlevel\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mpackage\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mlevel\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    127\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    128\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/_bootstrap.py\u001b[0m in \u001b[0;36m_gcd_import\u001b[0;34m(name, package, level)\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/_bootstrap.py\u001b[0m in \u001b[0;36m_find_and_load\u001b[0;34m(name, import_)\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/_bootstrap.py\u001b[0m in \u001b[0;36m_find_and_load_unlocked\u001b[0;34m(name, import_)\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/_bootstrap.py\u001b[0m in \u001b[0;36m_find_spec\u001b[0;34m(name, path, target)\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/_bootstrap_external.py\u001b[0m in \u001b[0;36mfind_spec\u001b[0;34m(cls, fullname, path, target)\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/_bootstrap_external.py\u001b[0m in \u001b[0;36m_get_spec\u001b[0;34m(cls, fullname, path, target)\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/_bootstrap_external.py\u001b[0m in \u001b[0;36mfind_spec\u001b[0;34m(self, fullname, target)\u001b[0m\n",
      "\u001b[0;32m/usr/lib/python3.10/importlib/_bootstrap_external.py\u001b[0m in \u001b[0;36m_path_stat\u001b[0;34m(path)\u001b[0m\n",
      "\u001b[0;31mKeyboardInterrupt\u001b[0m: "
     ]
    }
   ],
   "source": [
    "startdate = '2018-01-01'\n",
    "enddate = '2018-12-31'\n",
    "\n",
    "merge_all_datasets(startdate, enddate, 'merged_dataset', 'merged')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "b53fd131",
   "metadata": {},
   "source": [
    "##### 2018 ######"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "1899ef79",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "5th percentile of rlds (longwaveradiation) :  278.29998779296875\n",
      "95th percentile of rlds (longwaveradiation) :  400.20001220703125\n"
     ]
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 504x360 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "#2018 longwave\n",
    "ds_2018 = xr.open_dataset('merged_dataset/' + 'merged_20180101_20181231.nc')  \n",
    "\n",
    "\n",
    "ax = ds_2018['rlds'].plot.hist(figsize=(7,5), bins = 100)\n",
    "#ax.set_xlabel(\"length of the gaps [min]\")\n",
    "#ax.set_title('Histogram of gaps in the year 2018')\n",
    "\n",
    "print(\"5th percentile of rlds (longwaveradiation) : \",\n",
    "       np.percentile(ds_2018.rlds.values, 5))\n",
    "print(\"95th percentile of rlds (longwaveradiation) : \",\n",
    "       np.percentile(ds_2018.rlds.values, 95))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "c82a84ca",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "5th percentile of rlds (shortwaveradiation):  0.5\n",
      "95th percentile of rlds (shortwaveradiation):  863.0\n"
     ]
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 504x360 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "#2018 shortwaveradiation\n",
    "ds_2018 = xr.open_dataset('merged_dataset/' + 'merged_20180101_20181231.nc')  \n",
    "ds_2018['rlds_nan_flag'] = np.isnan(ds_2018.rlds)\n",
    "ds_2018 = ds_2018.where(ds_2018['rlds_nan_flag'] == False).dropna(how='all', dim='time')\n",
    "\n",
    "\n",
    "ax = ds_2018['rsds'].plot.hist(figsize=(7,5), bins = 100)\n",
    "ax.set_ylabel(\"frequency\")\n",
    "\n",
    "\n",
    "print(\"5th percentile of rlds (shortwaveradiation): \",\n",
    "       np.nanpercentile(ds_2018.rsds.values, 5))\n",
    "print(\"95th percentile of rlds (shortwaveradiation): \",\n",
    "       np.nanpercentile(ds_2018.rsds.values, 95))"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "7ee6e406",
   "metadata": {},
   "source": [
    "##### 2019 #####"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 91,
   "id": "5a6a89b9",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "5th percentile of rlds (longwaveradiation) :  261.0\n",
      "95th percentile of rlds (longwaveradiation) :  395.29998779296875\n"
     ]
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 504x360 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "#2019 longwave\n",
    "ds_2019 = xr.open_dataset('merged_dataset/' + 'merged_20190101_20191231.nc')  \n",
    "\n",
    "\n",
    "ax = ds_2019['rlds'].plot.hist(figsize=(7,5), bins = 100)\n",
    "ax.set_ylabel(\"frequency\")\n",
    "\n",
    "\n",
    "print(\"5th percentile of rlds (longwaveradiation) : \",\n",
    "       np.percentile(ds_2019.rlds.values, 5))\n",
    "print(\"95th percentile of rlds (longwaveradiation) : \",\n",
    "       np.percentile(ds_2019.rlds.values, 95))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "7465fd8f",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "5th percentile of rlds (shortwaveradiation):  0.10000000149011612\n",
      "95th percentile of rlds (shortwaveradiation):  831.2999877929688\n"
     ]
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 504x360 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "#2019 shortwaveradiation\n",
    "ds_2019 = xr.open_dataset('merged_dataset/' + 'merged_20190101_20191231.nc')  \n",
    "\n",
    "\n",
    "ax = ds_2019['rsds'].plot.hist(figsize=(7,5), bins = 100)\n",
    "#ax.set_xlabel(\"length of the gaps [min]\")\n",
    "#ax.set_title('Histogram of gaps in the year 2018')\n",
    "\n",
    "print(\"5th percentile of rlds (shortwaveradiation): \",\n",
    "       np.nanpercentile(ds_2019.rsds.values, 3))\n",
    "print(\"95th percentile of rlds (shortwaveradiation): \",\n",
    "       np.nanpercentile(ds_2019.rsds.values, 95))"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "b8fe3e2e",
   "metadata": {},
   "source": [
    "### Plot of the length of gaps of the data in a histogram ###\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 28,
   "id": "1c962c8b",
   "metadata": {
    "collapsed": true
   },
   "outputs": [
    {
     "ename": "KeyboardInterrupt",
     "evalue": "",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mKeyboardInterrupt\u001b[0m                         Traceback (most recent call last)",
      "\u001b[0;32m/tmp/ipykernel_1993196/124275217.py\u001b[0m in \u001b[0;36m<module>\u001b[0;34m\u001b[0m\n\u001b[1;32m      2\u001b[0m \u001b[0mstartdate\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m'2019-01-01'\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m      3\u001b[0m \u001b[0menddate\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m'2019-12-31'\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m----> 4\u001b[0;31m \u001b[0mdata_quality\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mstartdate\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0menddate\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'full_datasets'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'all_data'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'data_quality'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'gaps'\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m",
      "\u001b[0;32m/tmp/ipykernel_1993196/288377720.py\u001b[0m in \u001b[0;36mdata_quality\u001b[0;34m(startdate, enddate, inpath, inname, outpath, outname)\u001b[0m\n\u001b[1;32m    255\u001b[0m         \u001b[0mds\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0;34m'rsds_unrealistic_flag'\u001b[0m\u001b[0;34m]\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m(\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0;34m'time'\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mnp\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mfull\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mlen\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mrsds\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;32mFalse\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    256\u001b[0m         \u001b[0;32mfor\u001b[0m \u001b[0mi\u001b[0m \u001b[0;32min\u001b[0m \u001b[0mrange\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;36m0\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mlen\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mrsds\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 257\u001b[0;31m             \u001b[0;32mif\u001b[0m \u001b[0mds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mrsds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mvalues\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mi\u001b[0m\u001b[0;34m]\u001b[0m \u001b[0;34m>\u001b[0m \u001b[0;36m1500\u001b[0m \u001b[0;32mor\u001b[0m \u001b[0mds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mrsds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mvalues\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mi\u001b[0m\u001b[0;34m]\u001b[0m \u001b[0;34m<\u001b[0m \u001b[0;34m-\u001b[0m\u001b[0;36m5\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    258\u001b[0m                 \u001b[0mds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mrsds_unrealistic_flag\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mi\u001b[0m\u001b[0;34m]\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;32mTrue\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    259\u001b[0m                 \u001b[0mprint\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m'unrealistic rsds-values detected'\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mday\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mprint\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mrsds\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mvalues\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mi\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/dataarray.py\u001b[0m in \u001b[0;36mvalues\u001b[0;34m(self)\u001b[0m\n\u001b[1;32m    727\u001b[0m         \u001b[0mtype\u001b[0m \u001b[0mdoes\u001b[0m \u001b[0;32mnot\u001b[0m \u001b[0msupport\u001b[0m \u001b[0mcoercion\u001b[0m \u001b[0mlike\u001b[0m \u001b[0mthis\u001b[0m \u001b[0;34m(\u001b[0m\u001b[0me\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mg\u001b[0m\u001b[0;34m.\u001b[0m \u001b[0mcupy\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    728\u001b[0m         \"\"\"\n\u001b[0;32m--> 729\u001b[0;31m         \u001b[0;32mreturn\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mvariable\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mvalues\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    730\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    731\u001b[0m     \u001b[0;34m@\u001b[0m\u001b[0mvalues\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0msetter\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/variable.py\u001b[0m in \u001b[0;36mvalues\u001b[0;34m(self)\u001b[0m\n\u001b[1;32m    606\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0mvalues\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    607\u001b[0m         \u001b[0;34m\"\"\"The variable's data as a numpy.ndarray\"\"\"\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 608\u001b[0;31m         \u001b[0;32mreturn\u001b[0m \u001b[0m_as_array_or_item\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_data\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    609\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    610\u001b[0m     \u001b[0;34m@\u001b[0m\u001b[0mvalues\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0msetter\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/variable.py\u001b[0m in \u001b[0;36m_as_array_or_item\u001b[0;34m(data)\u001b[0m\n\u001b[1;32m    312\u001b[0m     \u001b[0mTODO\u001b[0m\u001b[0;34m:\u001b[0m \u001b[0mremove\u001b[0m \u001b[0mthis\u001b[0m \u001b[0;34m(\u001b[0m\u001b[0mreplace\u001b[0m \u001b[0;32mwith\u001b[0m \u001b[0mnp\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0masarray\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0monce\u001b[0m \u001b[0mthese\u001b[0m \u001b[0missues\u001b[0m \u001b[0mare\u001b[0m \u001b[0mfixed\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    313\u001b[0m     \"\"\"\n\u001b[0;32m--> 314\u001b[0;31m     \u001b[0mdata\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mnp\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0masarray\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mdata\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    315\u001b[0m     \u001b[0;32mif\u001b[0m \u001b[0mdata\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mndim\u001b[0m \u001b[0;34m==\u001b[0m \u001b[0;36m0\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    316\u001b[0m         \u001b[0;32mif\u001b[0m \u001b[0mdata\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mdtype\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mkind\u001b[0m \u001b[0;34m==\u001b[0m \u001b[0;34m\"M\"\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/indexing.py\u001b[0m in \u001b[0;36m__array__\u001b[0;34m(self, dtype)\u001b[0m\n\u001b[1;32m    652\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0m__array__\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mdtype\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0;32mNone\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    653\u001b[0m         \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_ensure_cached\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 654\u001b[0;31m         \u001b[0;32mreturn\u001b[0m \u001b[0mnp\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0masarray\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0marray\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mdtype\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mdtype\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    655\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    656\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0m__getitem__\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mkey\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/indexing.py\u001b[0m in \u001b[0;36m__array__\u001b[0;34m(self, dtype)\u001b[0m\n\u001b[1;32m    444\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0m__array__\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mdtype\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0;32mNone\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    445\u001b[0m         \u001b[0mkey\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mBasicIndexer\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mslice\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;32mNone\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0;34m)\u001b[0m \u001b[0;34m*\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mndim\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 446\u001b[0;31m         \u001b[0;32mreturn\u001b[0m \u001b[0mnp\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0masarray\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mkey\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mdtype\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0mdtype\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    447\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    448\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m~/.local/lib/python3.10/site-packages/xarray/core/indexing.py\u001b[0m in \u001b[0;36m__getitem__\u001b[0;34m(self, key)\u001b[0m\n\u001b[1;32m   1261\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   1262\u001b[0m     \u001b[0;32mdef\u001b[0m \u001b[0m__getitem__\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mself\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mkey\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m:\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m-> 1263\u001b[0;31m         \u001b[0marray\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mkey\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mself\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0m_indexing_array_and_key\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mkey\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m   1264\u001b[0m         \u001b[0;32mreturn\u001b[0m \u001b[0marray\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mkey\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m   1265\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;31mKeyboardInterrupt\u001b[0m: "
     ]
    }
   ],
   "source": [
    "#help(data_quality)\n",
    "startdate = '2019-01-01'\n",
    "enddate = '2019-12-31'\n",
    "data_quality(startdate, enddate, 'full_datasets', 'all_data', 'data_quality', 'gaps')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 29,
   "id": "ad98c2dc",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "Text(0.5, 1.0, 'Histogram of the length of gaps in the year 2018')"
      ]
     },
     "execution_count": 29,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 720x360 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "#2018\n",
    "df_2018 = pd.read_pickle('data_quality/' + 'gaps_20180101_20181231.pkl')  \n",
    "df_2018['length_of_gap']\n",
    "\n",
    "ax = df_2018['length_of_gap'].astype('timedelta64[m]').plot.hist(figsize=(10,5),bins = 60 )\n",
    "ax.set_xlabel(\"length of the length of gaps [min]\")\n",
    "ax.set_title('Histogram of the length of gaps in the year 2018')\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 30,
   "id": "96eed3d6",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "Text(0.5, 1.0, 'Histogram of the length of gaps in the year 2019')"
      ]
     },
     "execution_count": 30,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 720x360 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "#2019 \n",
    "df_2019 = pd.read_pickle('data_quality/' + 'gaps_20190101_20191231.pkl')  \n",
    "df_2019['length_of_gap']\n",
    "ax = df_2019['length_of_gap'].astype('timedelta64[m]').plot.hist(figsize=(10,5),bins = 60)\n",
    "ax.set_xlabel(\"length of the length of gaps [min]\")\n",
    "ax.set_title('Histogram of the length of gaps in the year 2019')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "78e44098",
   "metadata": {},
   "source": []
  },
  {
   "cell_type": "markdown",
   "id": "221c3343",
   "metadata": {},
   "source": [
    "### monthly timeseries with gaps per day ###"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 27,
   "id": "7de31c83",
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "/tmp/ipykernel_1993196/2143357989.py:28: SettingWithCopyWarning: \n",
      "A value is trying to be set on a copy of a slice from a DataFrame\n",
      "\n",
      "See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy\n",
      "  number_of_gaps[d] = number_of_gaps[d] + 1\n"
     ]
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 1080x360 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "gaps_per_day_in_a_year('2018')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 23,
   "id": "ac8f64f0",
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "/tmp/ipykernel_1993196/3510068300.py:28: SettingWithCopyWarning: \n",
      "A value is trying to be set on a copy of a slice from a DataFrame\n",
      "\n",
      "See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy\n",
      "  number_of_gaps[d] = number_of_gaps[d] + 1\n"
     ]
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 1080x360 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "gaps_per_day_in_a_year('2019')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "13efab58",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.10.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
