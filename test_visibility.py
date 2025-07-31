import numpy as np
np.float_ = np.float64

from wrf import getvar, to_np
from netCDF4 import Dataset

import metpy.calc as mpcalc
from metpy.units import units
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
from matplotlib.animation import FuncAnimation, PillowWriter

from NEW_visibility import calc_visibility
#from new_vis import calculate_visibility

def lowest_var(file, var, unit=False, surface=False, wind=False):
    data = getvar(file, var, units=unit) if unit else getvar(file, var, meta=False)
    if wind:
        data = data[0]
    data = to_np(data)

    return data[125:185,155:215] if surface else data[:,125:185,155:215]

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
img = ax.pcolormesh(lon, lat, np.zeros_like(lat), cmap='viridis', vmin=0, vmax=15)
cbar = plt.colorbar(img, ax=ax)
cbar.set_label('Visibility (mi)')
timestamp = ax.set_title(f'init {init_time} - {ens_num}')

vis_list = []
def update(frame):
    print(frame)
    time = wrf_files[frame].split('/')[-1].split('_')[-1]
    ds = Dataset(wrf_files[frame])
    pres_3d = lowest_var(ds, 'pres', unit='Pa')
#    temp_3d = lowest_var(ds, 'tk')
    temp_3d = lowest_var(ds, 'theta')
    u_3d = lowest_var(ds, 'ua')
    v_3d = lowest_var(ds, 'va')
    qv = lowest_var(ds, 'QVAPOR')
    qc = lowest_var(ds, 'QCLOUD')
    qr = lowest_var(ds, 'QRAIN')
    qi = lowest_var(ds, 'QICE')
    qs = lowest_var(ds, 'QSNOW')
    qg = lowest_var(ds, 'QGRAUP')
    czen = lowest_var(ds, 'COSZEN', surface=True)

    vis_grid = calc_visibility(
            pres_3d, temp_3d, u_3d, v_3d, qv, qc,
            qr, qi, qs, qg, czen)
 
#    vis_grid = calculate_visibility(qv, qc, qr, qi, qs, temp_3d, pres_3d)
    vis_list.append(vis_grid)

    img.set_array(vis_grid.ravel())
    timestamp.set_text(f'{time}')
    ds.close()

    return [img, timestamp]

ani = FuncAnimation(fig, update, frames=len(wrf_files), blit=False)
ani.save('vis_anim.gif', writer=PillowWriter(fps=5))

vis_all = np.stack(vis_list, axis=0)
np.save(f'{ens_num}_vis.npy', vis_all) #.filled(np.nan))
print(np.nanmin(vis_list))
