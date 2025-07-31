import numpy as np
np.float_ = np.float64

from wrf import getvar, to_np
from netCDF4 import Dataset

import os
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import quantities
import time
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.animation import FuncAnimation, PillowWriter

from SSP import snow_squall_param

def lowest_var(file, var, unit=False, surface=False, wind=False):
    data = getvar(file, var, units=unit) if unit else getvar(file, var, meta=False)
    if wind:
        data = data[0]
    data = to_np(data)

    return data[125:185,155:215] if surface else data[:,125:185,155:215]

def get_windspeed():
    return None

path = '/scratch4/STI/coastal/Carolyn.Emerson/lapenta_ww/wofs_WY_snow/0000/ENS_1/'

init_time = path.split('/')[6]
ens_num = path.split('/')[7]

wrf_files = sorted([os.path.join(path, f) for f in os.listdir(path)])

# first file
ds = Dataset(wrf_files[0], mode='r')
lat = lowest_var(ds, 'lat', surface=True)
lon = lowest_var(ds, 'lon', surface=True)

# set up plot
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ax.coastlines()
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5)
img = ax.pcolormesh(lon, lat, np.zeros_like(lat), cmap='viridis', vmin=0, vmax=6)
cbar = plt.colorbar(img, ax=ax)
cbar.set_label('Snow Squall Parameter (unitless)')
timestamp = ax.set_title(f'init {init_time} - {ens_num}')

snsq_list = []
def update(frame):
    time = wrf_files[frame].split('/')[-1].split('_')[-1]
    ds = Dataset(wrf_files[frame])
    pres_3d = lowest_var(ds, 'pres', unit='Pa')
    temp_3d = lowest_var(ds, 'tk')
    height_3d = lowest_var(ds, 'height_agl', unit='m')
    wind_3d = lowest_var(ds, 'uvmet_wspd_wdir', wind=True)
    td_3d = lowest_var(ds, 'td', unit='K')
    rh_3d = lowest_var(ds, 'rh')

    pres_sfc = lowest_var(ds, 'PSFC', surface=True)
    temp_sfc = lowest_var(ds, 'T2', surface=True)
    wind_sfc = lowest_var(ds, 'uvmet10_wspd_wdir', surface=True, wind=True)
    td_sfc = lowest_var(ds, 'td2', unit='K', surface=True)
    rh_sfc = lowest_var(ds, 'rh2', surface=True)
    msl_sfc = lowest_var(ds, 'ter', surface=True)

    print(f'{frame}')
    snsq_grid = snow_squall_param(
            pres_3d, temp_3d, td_3d, rh_3d, height_3d, wind_3d,
            msl_sfc, pres_sfc, temp_sfc, td_sfc, rh_sfc, wind_sfc)
    snsq_list.append(snsq_grid)

#    print(snsq_grid)
#    print(np.nanmean(snsq_grid))
    img.set_array(snsq_grid.ravel())
    timestamp.set_text(f'{time}')
    ds.close()

    return [img, timestamp]

ani = FuncAnimation(fig, update, frames=len(wrf_files), blit=False)
ani.save('SSP_anim.gif', writer=PillowWriter(fps=5))

snsq_all = np.stack(snsq_list, axis=0)
np.save(f'{ens_num}_snsq.npy', snsq_all)


print('done')
