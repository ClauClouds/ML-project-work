################################################################################
# this script was written by Jose, I made slight adjustments, but mainly performance wise. This script contains all the functions used to calculate Ka-X-spectra and noise, aswell as the edges of the Ka-spectra
###############################################################################

import os
import numpy as np
import pandas as pd
import xarray as xr
import glob as gb
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def getCalSpect(dataSet, ch):
    '''
    author: Jose dias Nieto, readapted by Claudia Acquistapace
    date: december 2022
    goal: function to calculate the calibrated doppler spectra
    input: 
        dataSet = dataset xarray containing Doppler spectra extracted from ncdf files 
        ch = frequency of the radar 
    output:
        calSPC = calibrated Doppler spectra
    '''
# function from Jose, I added notes to that which might not be correct...    
#ch = o or x (co channel or cross channel)
    
    rangeLen = dataSet.range.shape[0]
    try:
        timeLen = dataSet.time.shape[0]
    except:
        timeLen = 1
    
    rangeData = dataSet.range.values.reshape(1, rangeLen, 1)
    #-- radar Constant related to 5km Height and 200ns pulses (description from nc-file)
    radConst = dataSet.RadarConst.values.reshape(timeLen,1,1)
    #-- SNRCorFa: Signal-to-noise-factor corrected based on Receiver calibration
    SNRCorFa = dataSet['SNRCorFaC'+ch].values.reshape(timeLen,rangeLen,1)
    #-- Doppler Spectrum
    SPC = dataSet['SPCc'+ch].values
    #-- Noise power
    npw = dataSet.npw1.values.reshape(timeLen,1,1)
    #-- provided by Stefan    
    calSPC = radConst*SNRCorFa*(rangeData**2/5000.**2)*SPC/npw
    
    return calSPC


def getCalNoise(dataSet, ch):
    '''
    author: Jose dias Nieto, readapted by Claudia Acquistapace
    date: december 2022
    goal: function to calculate the noise level of the Doppler spectra
    input: 
        dataSet = dataset xarray containing Doppler spectra extracted from ncdf files 
        ch = frequency of the radar 
    output:
        calNoise = noise level
    '''    
    rangeLen = dataSet.range.shape[0]
    try:
        timeLen = dataSet.time.shape[0]
    except:
        timeLen = 1
    
    rangeData = dataSet.range.values.reshape(1, rangeLen, 1)
    radConst = dataSet.RadarConst.values.reshape(timeLen,1,1)
    SNRCorFa = dataSet['SNRCorFaC'+ch].values.reshape(timeLen,rangeLen,1)
    #-- ???? What does HSD mean????
    HSD = dataSet['HSDc'+ch].values.reshape(timeLen, rangeLen, 1)
    npw = dataSet.npw1.values.reshape(timeLen,1,1)
    
    calNoise = radConst*SNRCorFa*(rangeData**2/5000.**2)*HSD/npw
    
    return calNoise


def reorderSpec(dopplerVal, specData, dpMin, dpMax):
    '''
    author: Jose dias Nieto, readapted by Claudia Acquistapace
    date: december 2022
    goal: function to reorder Doppler velocity array that is not provided in growing order
    input: 
        dopplerVal = array containing values of Doppler velocity
        specData = Doppler spectra
        dpMin = max Doppler velocity
        dpMax = min Doppler velocity
        ch = frequency of the radar 
    output:
        newDoppler = new growing Doppler velocity array
        newSpec = reordered Doppler spectrum
    '''    

    # Joses skript!    
    #-- reorder Spectrum around maximum value (-> switch left and right side of spectrum)
    dpIndex = np.arange(len(dopplerVal))
    
    maxDpIndex = dpIndex[dopplerVal >= 0].max()
    zeroIndex = dpIndex[dopplerVal >= 0].min()
    
    newDoppler = np.ones_like(dopplerVal)
    newSpec = np.ones_like(specData) * np.nan
    
    newDoppler[len(newDoppler) - maxDpIndex-1:] = dopplerVal[:maxDpIndex+1]
    newDoppler[:len(newDoppler) - maxDpIndex-1] = dopplerVal[maxDpIndex+1:]
    
    newSpec[:,:,len(newDoppler) - maxDpIndex-1:] = specData[:,:,:maxDpIndex+1]
    newSpec[:,:,:len(newDoppler) - maxDpIndex-1] = specData[:,:,maxDpIndex+1:]

    newDoppler = newDoppler*(-1)
    try:
        doppLimIndMin = dpIndex[newDoppler>=dpMax].max()
    except:
        doppLimIndMin = 0
    try:
        doppLimIndMax = dpIndex[newDoppler<=dpMin].min()
    except:
        doppLimIndMax = len(dopplerVal)
    
    newDoppler = newDoppler[doppLimIndMin:doppLimIndMax+1]
    newSpec = newSpec[:,:,doppLimIndMin:doppLimIndMax+1]
   
    return newDoppler, newSpec


def getMaxAndMinSpecVel(spectra, velArr):
    '''
    author: Jose dias Nieto, readapted by Claudia Acquistapace
    date: december 2022
    goal: function to calculate max and min velocity
    input: 
        spectra = 
        velArr = array containing values of Doppler velocity
    output:
        maxVel = max velocity
        minVel = min velocity
    ''' 
    # Joses script!
    #-- calculates max and min velocity    
    velMatrix = (spectra/spectra)*velArr
    maxVel = np.nanmax(velMatrix, axis=1)
    minVel = np.nanmin(velMatrix, axis=1)

    return maxVel, minVel


def loadKaData(dataPathKa, wTime):
    '''
    author: Jose dias Nieto, readapted by Claudia Acquistapace
    date: december 2022
    goal: function to r4ead ka band datasets
    input: 
        dataPathKa = string containing path to the data of Ka band radar
        wTime = time array
    output:
        kaData = ka band radar dataset
    ''' 
#- Joses script!
#- opens dataset from Ka band at that specific hour
    if pd.to_datetime(wTime).hour % 2 == 0:

        fileId = pd.to_datetime(wTime).strftime('%Y%m%d_%H')
    
    else:
    
        fileId = pd.to_datetime(wTime) - pd.to_timedelta('1H')
        fileId = fileId.strftime('%Y%m%d_%H')
    
    filePathKa = os.path.join(dataPathKa, fileId+'*.znc')
    fileListW = gb.glob(filePathKa)

    kaData = xr.open_dataset(fileListW[0])
    kaData.time.values = kaData.time.values + kaData.microsec.values*10**(-6)
    kaData.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC'
    kaData = xr.decode_cf(kaData)
    
    return kaData


def removSuspData(calNoise, newSpecKa, noiseCutoff):
        '''
    author: Jose dias Nieto, readapted by Claudia Acquistapace
    date: december 2022
    goal: function to remove suspicious spectra data, using edge velocities at 3dB above noise?
    input: 
        calNoise = calibrated noise level
        newSpecKa = calibrated Doppler spectra
        noiseCutoff = noise cut off level
    output:
        newSpecKa = new spectra 
    ''' 
#- joses script, use edge velocities at 3db above noise?

    noiseLevelArr = 10*np.log10(calNoise[0,:,0])+noiseCutoff

    for index, noiseLevel in enumerate(noiseLevelArr):

        specRow = newSpecKa[0][index]
        specRow[specRow < noiseLevel] = np.nan
        newSpecKa[0][index] = specRow
        
    return newSpecKa


def getMaskProf(kaData, newSpecKa, newDopplerKa, minHeight, Nyquist):
    '''
    author: Jose dias Nieto, readapted by Claudia Acquistapace
    date: december 2022
    goal: function to remove suspicious spectra data, using edge velocities at 3dB above noise?
    input: 
        kaData, newSpecKa, newDopplerKa, minHeight, Nyquist
    output:
        newSpecKa = new spectra 
    ''' 
#- Joses script!
#- calculate max and min velocities    
    maxVelKa, minVelKa = getMaxAndMinSpecVel(newSpecKa[0], newDopplerKa)

    zeros = np.zeros_like(kaData.range.values)
    maxVelKa = np.nansum(np.array([zeros,maxVelKa]), axis=0)
    minVelKa = np.nansum(np.array([zeros,minVelKa]), axis=0)

    maxVelKa[kaData.range.values<minHeight] = np.nan
    minVelKa[kaData.range.values<minHeight] = np.nan
    maxVelKa[maxVelKa>= Nyquist] = Nyquist + 10
    minVelKa[minVelKa<=-Nyquist] = Nyquist - 10

    return maxVelKa, minVelKa



def f_calc_spec_ka(spec_data, channel):
    '''function to calculate calibrated Doppler spectra, reorder its bins and convert it in dB units
    input: spectra xarray dataset
           channel (options 'o', 'x')
    output: 
    '''
    # now most of the functions used here are directly used from Jose
    #for channel in ['o','x']:
    
    #calculate calibrated spectra
    calSpect = kaLib.getCalSpect(spec_data, channel)
    
    # calculate calibrated noise
    calNoise = kaLib.getCalNoise(spec_data, channel)
    
    # calculate spec - noise
    specCalNoiCorrKa = calSpect - (np.ones_like(calSpect)*calNoise)
    #-- reorder the spectrum because the ka spectrum is stored as 0,..,10,-10,..0 
    newDopplerKa, newSpecKa = kaLib.reorderSpec(spec_data.doppler.values,
                                                specCalNoiCorrKa,
                                                -10.56824, 
                                                10.56824)# NyquistVelocity

    # calculation of log Ze
    newSpecKa = 10*np.log10(newSpecKa)
    
    #removing suspicious data
    newSpecKa = kaLib.removSuspData(calNoise, newSpecKa, 3)
    
    SpecKa = xr.DataArray(newSpecKa, 
                          dims=('time','range','doppler'),
                          coords={'time':spec_data.time.values,'range':spec_data.range, 'doppler':newDopplerKa})
    return(SpecKa, newSpecKa)


# set the colormap and centre the colorbar
class MidpointNormalize(mpl.colors.Normalize):
    """Normalise the colorbar."""
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    

def f_extract_date(date):
    '''reads date in the format yyyymmdd and returns yyyy, mm, dd'''
    yyyy = date[0:4]
    mm = date[4:6]
    dd = date[6:8]
    
    return(yyyy, mm, dd)

def f_read_file_list(yyyy, mm, dd, hour):
    
    '''reads the file list associated to the selected date and returns it'''
    # list all files ncdf for the day
    file_list = glob.glob('/data/obs/site/jue/joyrad35/'+yyyy+'/'+mm+'/'+dd+'/'+yyyy+mm+dd+'_'+hour+'*.znc')
    
    #print(file_list)
    # remove files with ppi.znc extensions
    list_ppi = np.sort(glob.glob('/data/obs/site/jue/joyrad35/'+yyyy+'/'+mm+'/'+dd+'/'+yyyy+mm+dd+'_'+hour+'*ppi.znc'))

    # remove ppi files from main list 
    for i_list,el in enumerate(list_ppi):
        file_list.remove(el)

    # sorting files remained in the final list
    file_list = np.sort(file_list)
    print(file_list)
    return(file_list)

def f_plot_Doppler_spectrogram(spec_slice_10_sec, doppler, range_height, date, hour, minute, path_out):
    '''Function to plot the height spectrogram mean over 10 minutes around the satellite time stamp'''
    
        # plot settings
    dict_plot_settings = {
        'plot_ticks'   :16,
        'labelsizeaxes':16,
        'fontSizeTitle':16,
        'fontSizeX'    :16,
        'fontSizeY'    :16,
        'cbarAspect'   :16,
        'fontSizeCbar' :16,
        'rcparams_font':['Tahoma'],
        'savefig_dpi'  :100,
        'font_size'    :16,
        'grid'         :True}


    # plots settings defined by user at the top
    labelsizeaxes   = dict_plot_settings['labelsizeaxes']
    fontSizeTitle   = dict_plot_settings['fontSizeTitle']
    fontSizeX       = dict_plot_settings['fontSizeX']
    fontSizeY       = dict_plot_settings['fontSizeY']
    cbarAspect      = dict_plot_settings['cbarAspect']
    fontSizeCbar    = dict_plot_settings['fontSizeCbar']
    rcParams['font.sans-serif'] = dict_plot_settings['rcparams_font']
    matplotlib.rcParams['savefig.dpi'] = dict_plot_settings['savefig_dpi']
    plt.rcParams.update({'font.size':dict_plot_settings['font_size']})
    grid = dict_plot_settings['grid']
    #matplotlib.rc('xtick', labelsize=dict_plot_settings['plot_ticks'])  # sets dimension of ticks in the plots
    #matplotlib.rc('ytick', labelsize=dict_plot_settings['plot_ticks'])  # sets dimension of ticks in the plots
    plt.rcParams["figure.figsize"]=(20, 16)

    
    # plot figure for ML input
    ymin = 100.
    ymax = 12000.
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,10))
    rcParams['font.sans-serif'] = ['Tahoma']
    matplotlib.rcParams['savefig.dpi'] = 100
    plt.gcf().subplots_adjust(bottom=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    mesh1 = ax.pcolormesh(doppler, range_height, spec_slice_10_sec.data,\
                                    cmap='jet', vmin=-80., vmax=10., rasterized=True)
    #cbar = fig.colorbar(mesh1, ax=ax)
    #cbar.set_label(label='Power [dB]',  size=20)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-10, 10.)
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_major_locator(ticker.NullLocator())
    #ax.axvline(x=0., linestyle='--', color='grey')
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    #ax.spines["bottom"].set_linewidth(3)
    #ax.spines["left"].set_linewidth(3)
    #ax.set_xlabel('Doppler velocity [ms$^{-1}$]')
    #ax.set_ylabel('Height [m]')

    #ax_dict['D'].axhline(y=radar_slice.cloud_base.values, linestyle='--', color='grey')
    #ax.set_title('10 min mean spectrogram', fontweight='black',\
    #             fontsize=fontSizeX, transform=ax.transAxes)
    fig.savefig(path_out+date+'_'+hour+minute+'_DS.png')
    
    return(1)


def f_save_to_ncdf(spec_slice_10_sec, doppler, range_height, date, hour, minute, path_out):
    '''Function to save ncdf files of doppler spectra for further postprocessing 
    ( calculation of mean, std and normalization of images for input to machine learning'''
    
    
    # constructing time stamp to associate to the data
    yy, mm, dd = f_extract_date(date)
    time_value = datetime(int(yy),int(mm),int(dd),int(hour), int(minute))
    print(time_value)

     # saving data in CF convention dataset
    Ka_radar_data = xr.Dataset(
        data_vars={
            'spectrum_linear':(('height','doppler'), (10**(spec_slice_10_sec/10)), \
                               {'long_name': 'Doppler spectrum', 'units':'mm6 m-3'}),
            'spectrum_db':(('height','doppler'), spec_slice_10_sec, \
                           {'long_name': 'Doppler spectrum', 'units':'dB'}),            
        },
        coords={
            "time": (('time',), [time_value], {"axis": "T","standard_name": "time"}), # leave units intentionally blank, to be defined in the encoding
            "height": (('height',), range_height, {"axis": "Z","positive": "up","units": "m", "long_name":'radar_range_height'}),
            'doppler': (('doppler',), doppler, { "axis": "D", "standard_name": "doppler velocity"}),

        },
        attrs={'CREATED_BY'     : 'Claudia Acquistapace',
            'CREATED_ON'       : str(datetime.now()),
            'FILL_VALUE'       : 'NaN',
            'PI_NAME'          : 'Claudia Acquistapace',
            'PI_AFFILIATION'   : 'University of Cologne (UNI), Germany',
            'PI_ADDRESS'       : 'Institute for geophysics and meteorology, Pohligstrasse 3, 50969 Koeln',
            'PI_MAIL'          : 'cacquist@meteo.uni-koeln.de',
            'DATA_DESCRIPTION' : 'Doppler spectra profiles at selected time stamps',
            'DATA_DISCIPLINE'  : 'Atmospheric Physics - Remote Sensing Radar Profiler',
            'DATA_GROUP'       : 'Experimental;Profile;site;Joyce',
            'DATA_LOCATION'    : 'Juelich, Germany',
            'DATA_SOURCE'      : 'Ka band Doppler spectra data postprocessed',
            'DATA_PROCESSING'  : 'https://github.com/ClauClouds/ML_project_Juelich/',
            'INSTRUMENT_MODEL' : '36 GHz (W-band) radar',
            'COMMENT'          : '' }
    )



    # assign istrument id
    instrument_id = xr.DataArray("Ka_band_radar",dims=(),attrs={"cf_role": "profiles"},)
    Ka_radar_data = Ka_radar_data.assign({"instrument": instrument_id,})

    # assign additional attributes following CF convention
    Ka_radar_data = Ka_radar_data.assign_attrs({
            "Conventions": "CF-1.8",
            "title": Ka_radar_data.attrs["DATA_DESCRIPTION"],
            "institution": 'University of Cologne (UNI), Germany',
            "history": "".join([
                "source: " + Ka_radar_data.attrs["DATA_SOURCE"] + "\n",
                "processing: " + Ka_radar_data.attrs["DATA_PROCESSING"] + "\n",
                " adapted to enhance CF compatibility\n",
            ]),  # the idea of this attribute is that each applied transformation is appended to create something like a log
            "featureType": "trajectoryProfile",
        })

    # storing ncdf data
    Ka_radar_data.to_netcdf(path_out+'/ncdf/'+date+hour+minute+'00.nc', \
                            encoding={"spectrum_linear":{"dtype": "f4", "zlib":True, "complevel":9},\
                                     "spectrum_db": {"dtype": "f4", "zlib": True, "complevel":9}, \
                                     "doppler": {"zlib": True, "complevel":9}, \
                                     "time": {"units": "seconds since 2020-01-01", "dtype": "i4"}})

    
