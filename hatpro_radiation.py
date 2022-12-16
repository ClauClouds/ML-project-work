from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import os, glob, sys, shutil
import math
import cv2
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

# path to radiation file and filename
loc = '/data/hatpro/jue/data/radiation/2013/06/'
file = 'CMPT_20130606.dat'
filename = loc+file

# reading variables from .dat file
timestamp,short_wave, error = [],[],[]
with open(filename) as input_data:
    temp = input_data.readlines()
    n = len(temp) #int(temp.pop(0))
    matrix = [x.strip().split(",") for x in temp[9:n]]
    for i in range(0,len(matrix)):
        timestamp.append(float(matrix[i][0]))
        short_wave.append(float(matrix[i][1]))
        error.append(float(matrix[i][2]))
    input_data.close()

# converting variables read in a pandas dataframe
df = pd.DataFrame()
df['timestamp'],df['short_wave'],df['error'] = '','',''
df['timestamp'],df['short_wave'],df['error'] = timestamp, short_wave, error
df['datetime'] = pd.to_datetime(df['timestamp'],unit='s', origin =pd.Timestamp('2000-01-01 00:00'))
df=df.set_index('datetime') #raw data available in every 1 sec
df1 = df.resample('5min').mean() # since sat data in 5 mins thats why resampling
