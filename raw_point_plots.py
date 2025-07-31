import numpy as np
np.float_ = np.float64

from wrf import getvar, to_np
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from NEW_visibility import calc_visibility

lat, lon = 41.25, -103.5

def lowest_var(file, var, points=None, unit=False, surface=False, wind=False):
    data = getvar(file, var, units=unit) if unit else getvar(file, var, meta=False)
    if wind:
        data = data[0]
    data = to_np(data)

    if surface:
        if points:
            return data[points[0]-1:points[0]+1, points[1]-1:points[1]+1]
        else:
            return data

    else:
        if points:
            return data[:,points[0]-1:points[0]+1, points[1]-1:points[1]+1]
        else:
            return data

# Path to your .nc files
directory = f'wofs_WY_snow/0000/'
cutoff_time = datetime.strptime('20200131_04:05:00', '%Y%m%d_%H:%M:%S')

ens_files = []
for foldername in os.listdir(directory): # loop through ensembles
    path = os.path.join(directory, foldername)
    time_files = sorted(os.listdir(path))

    time_paths = []
    for file in time_files:
        timestamp_str = file.split('_')[2] + '_' + file.split('_')[3]
        file_time = datetime.strptime(timestamp_str, '%Y-%m-%d_%H:%M:%S')

        if file_time<cutoff_time:
            time_paths.append(directory+'/'+foldername+'/'+file)
        
    time_files = np.array(time_paths)
    ens_files.append(time_files)

ens_files = np.array(ens_files)
n_ens, n_times = ens_files.shape


timestr = '202001310000'
time = datetime.strptime(timestr, '%Y%m%d%H%M')
timestamps = [time+timedelta(minutes=5*i) for i in range(n_times)]

fig, ax = plt.subplots(figsize=(9,7))

ens_means = []
for ens in [3,5,15,17]:
    time_vals = []
    for time in range(n_times):
        ds = Dataset(ens_files[ens, time], mode='r')
        lats = lowest_var(ds, 'lat', surface=True)
        lons = lowest_var(ds, 'lon', surface=True)

        # Compute squared distance
        dist = (lats - lat)**2 + (lons - lon)**2
        points = np.unravel_index(np.argmin(dist), dist.shape)

        temp = lowest_var(ds, 'T2', points, surface=True)

        ds.close()
        vals = temp[1, 1]
        
        vals = (vals-273.15)*(9/5) + 32
        time_vals.append(vals)
    ax.plot(timestamps, time_vals, linestyle='dotted', label=f'Ens {ens}')
    ens_means.append(time_vals)

ens_means = np.array(ens_means).mean(axis=0)
ax.plot(timestamps, ens_means, linestyle='solid', color='red', label='Ens. Mean.')
ax.set_xlim(timestamps[0], timestamps[-1])
ax.set_ylim(26.5,35)
ax.legend()

time1 = datetime(2020,1,31,2,0)
time2 = datetime(2020,1,31,3,0)

ax.axvline(time1, color='darkgray', linestyle='--', linewidth=0.8)
ax.axvline(time2, color='darkgray', linestyle='--', linewidth=0.8)

fig.text(0.13, 0.9, '4/18 WoFS Ens. 2m Temperature', fontsize=10, fontweight='bold',
        ha='left', va='top')
fig.text(0.9, 0.9, 'Init: 2020-01-31 0000 UTC', fontsize=10, ha='right', va='top')
ax.set_ylabel('2m Temperature (°F}')

formatter = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(formatter)

plt.savefig(f'tempF.png')

