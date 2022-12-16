#!/usr/bin/env python
# coding: utf-8
'''
# # Script example to read and plot doppler spectra from JOYCE supersite

# author: Claudia Acquistapace \\
# date: 28.11.2022
#
# general info: to select a case study please check the data browser:
# https://atmos.meteo.uni-koeln.de/~hatpro/dataBrowser/dataBrowser4.html
#
'''



import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import glob
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from matplotlib import rcParams
from warnings import warn
import datetime as dt
from scipy import interpolate
import custom_color_palette as ccp
import matplotlib as mpl
import os.path
import itertools
import os.path
import metpy.calc as mpcalc
from metpy.units import units
import metpy
import netCDF4 as nc4
import kaLib


# In[23]:


# define the date to process based on satellite data input
dates = ['20221128', '20221128']
hours = ['03', '04']
minutes = ['16', '10']

# set the number of time steps over which to average the spectra for deriving mean profile
'''time resolution of radar obs is 2 seconds, hence 300 time stamps == 5 minutes,
600 time stamps = 10 minutes. Consider that we add and subtract delta from the closest
time stamp, so we need to split the time interval in 2.'''
delta = 150


# path for joyce Doppler spectra data and path for output plot
path_out = '/net/ostro/ML_work_joyce/'

# counting number of images in input
Number_images = len(dates)

for ind_image, date in enumerate(dates):

    print('processing '+date)
    # reading corresponding hour and minute to the selected image
    hour = hours[ind_image]
    minute = minutes[ind_image]

    # estracting strings of year, month, day, hour, minute, second
    yyyy, mm, dd = f_extract_date(date)

    # defining time stamp to look for:
    df = pd.DataFrame({'year':[int(yyyy)],                        'month':[int(mm)],                        'day':[int(dd)],                        'hour':[int(hour)],                        'minute':[int(minute)]})
    TimeVal = pd.to_datetime(df)[0]


    # read file list of joyce data for the day
    file_list_spec_joyce = f_read_file_list(yyyy, mm, dd, hour)

    # read data in a spectra file
    spec_data = xr.open_mfdataset(file_list_spec_joyce)

    # reading elevation and azimut to select only zenith pointing data
    elevation = spec_data.elv
    azi = spec_data.azi

    # define a flag for filtering data
    filter_zenith = np.ones(len(spec_data.time.values))
    filter_zenith[np.where(elevation != 90.)[0]] = 0.

    # filtering data
    if len(spec_data.time.values[filter_zenith == 1]) == 0:
        print('skip date:'+ date+'_'+hour + minute)
    else:
        # select only zenith values
        spec_data_zenith = spec_data.sel(time=spec_data.time.values[filter_zenith == 1])

        # calculating calibrated spectra for co and cx channels
        new_spec_co_DA, specKa_co = f_calc_spec_ka(spec_data_zenith, 'o')
        #new_spec_cx_DA, specKa_cx = f_calc_spec_ka(spec_data_zenith, 'x')

        #- calculate radar reflectivity (Ze) and mean DOppler velocity (mdv) for Ka band radar
        #ZeKa = (10**(new_spec_co_DA/10)).sum(dim='doppler')
        #mdvKa = ((10**(new_spec_co_DA/10)*(new_spec_co_DA['doppler'])).sum(dim='doppler'))/ZeKa

        # read variables to plot and time array
        timeRadar = pd.to_datetime(new_spec_co_DA.time.values, unit='s')
        start_time = timeRadar[0]
        end_time = timeRadar[-1]

        for ind in range(len(timeRadar)):
            if timeRadar[ind] > TimeVal:
                print(ind)
                break

        print('selected index', ind)
        print('time in the time array', timeRadar[ind])
        print('value selected', TimeVal)

        # averaging spectra over 5 minutes
        spec_slice_10_sec = np.nanmedian(new_spec_co_DA[ind-delta:ind+delta, :, :], axis=0)

        #plotting the doppler spectrogram
        doppler = new_spec_co_DA.doppler
        range_height = new_spec_co_DA.range
        plotDone = f_plot_Doppler_spectrogram(spec_slice_10_sec, doppler, range_height, date, hour, minute, path_out)

        # saving variable Doppler spectra (dB and linear) and its coordinates in ncdf
        saved_data = f_save_to_ncdf(spec_slice_10_sec, doppler, range_height, date, hour, minute, path_out)
        print('data stored in ncdf')
