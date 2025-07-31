import numpy as np
np.float_ = np.float64

from wrf import getvar, to_np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import pandas as pd
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from datetime import datetime, timedelta
from plotting_subroutines import apply_max_or_min_filter, _gaussian_filter
from NEW_visibility import calc_visibility

refl_thresh = 25
max_size = 3
gauss_size = 1

directory = 'wofs_WY_snow/0000'
run_time = directory.split('/')[1]

roads_shp = 'roads/tl_2024_us_primaryroads.shp'
roads_feature = ShapelyFeature(
        Reader(roads_shp).geometries(),
        crs=ccrs.PlateCarree(),
        edgecolor='gray',
        facecolor='none')

def lowest_var(file, var, unit=False, surface=False, wind=False):
    data = getvar(file, var, units=unit) if unit else getvar(file, var, meta=False)
    if wind:
        data = data[0]
    data = to_np(data)

    return data[115:195,155:225] if surface else data[:,115:195,155:225]

def get_reflectivity(ds):
    refl = lowest_var(ds, 'dbz')[0,:,:]
    return refl

# get file paths
cutoff_time = datetime.strptime('20200131_04:05:00', '%Y%m%d_%H:%M:%S')

ens_files = []
for foldername in os.listdir(directory):
    path = os.path.join(directory, foldername)
    time_files = sorted(os.listdir(path))

    time_paths = []
    for file in time_files:
        timestamp_str = file.split('_')[2]+'_'+file.split('_')[3]
        file_time = datetime.strptime(timestamp_str, '%Y-%m-%d_%H:%M:%S')

        if file_time < cutoff_time:
            time_paths.append(directory+'/'+foldername+'/'+file)
    time_files = np.array(time_paths)

    ens_files.append(time_files)

ens_files = np.array(ens_files)
n_ens, n_times = ens_files.shape

# first file
ds = Dataset(ens_files[0,0], mode='r')
lats = lowest_var(ds, 'lat', surface=True)
lons = lowest_var(ds, 'lon', surface=True)
ny, nx = lats.shape
ds.close()

# create a disk memmap for probability fields
probs_map = np.lib.format.open_memmap(
        'probs.npy', mode='w+',
        dtype='float32', shape=(n_times, ny, nx)
        )

for time in range(n_times):
    print(f'Time: {time}')

    count = np.zeros((ny, nx), dtype=int)

    for ens in ens_files:
        ds = Dataset(ens[time], mode='r')
        var = lowest_var(ds, 'uvmet10_wspd_wdir', unit='mi h-1', surface=True, wind=True)

        ds.close()

        var_max = apply_max_or_min_filter(data=var, max_size=max_size)
        mask = (var_max >= refl_thresh)
        count += mask.astype(int)

    prob_raw = (count / float(n_ens)) * 100 # raw probability at each grid point
    probs_map[time, :, :] = _gaussian_filter(data=prob_raw, size=gauss_size)

probs_map.flush()
del probs_map

probs_map = np.lib.format.open_memmap(
        'probs.npy', mode='r'
        )

levels = np.arange(0,101,10)

# Light‐to‐dark blues for 0–40; light‐to‐dark reds for 50–100
colors = [
    "#d0e8ff", "#a1d1ff", "#72baff", "#439fff",  # 0–10,10–20,20–30,30–40
    "#ffe5e5", "#ffb2b2", "#ff7f7f", "#ff4c4c",
    "#ff1919", "#cc0000"                         # 40–50,...,90–100
]

cmap = ListedColormap(colors)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

artists = {}

fig, ax = plt.subplots(figsize=(12,9),subplot_kw={'projection':ccrs.PlateCarree()})

artists['img'] = ax.contourf(lons, lats, probs_map[0], levels=levels[1:], cmap=cmap, norm=norm,
        alpha=0.6, transform=ccrs.PlateCarree())
artists['cs'] = ax.contour(lons, lats, probs_map[0], levels=levels[1:], colors=[cmap(norm(l)) for l in levels],
        linewidths=0.5, transform=ccrs.PlateCarree())

cbar = plt.colorbar(artists['img'], ax=ax, pad=0.02, shrink=0.6)
cbar.set_label('Ensemble Probability', fontsize=9)

fig.text(0.13, 0.86, 'Ens. Prob. of 10m Wind Speed > 25 mph (3 km Neigh.)',
            fontsize=10, fontweight='bold', ha='left', va='top')

fig.text(0.13, 0.84, 'Ens. 2 Snow Squall Reflectivity (20 dBz, contour)', fontsize=10, fontweight='normal',
        ha='left', va='top')

init_str = 'Init: 2020-01-31 0000 UTC'
fig.text(0.77, 0.86, init_str, fontsize=8, fontweight='normal', ha='right', va='top')

ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='gray',linewidth=1)
ax.coastlines()
ax.set_extent([-104.5, -102.5, 40.5, 42])
ax.add_feature(roads_feature, linewidth=0.45, edgecolor='darkgray')

curr_str = 'Valid: 2020-01-31 0000 UTC'
valid_time = fig.text(0.77, 0.84, curr_str, fontsize=8, fontweight='normal', ha='right', va='top')

time = curr_str.split(' ')[1] + ' ' + curr_str.split(' ')[2]
time = datetime.strptime(time, '%Y-%m-%d %H%M')
timestamps = [time+timedelta(minutes=5*i) for i in range(n_times)]

fig.canvas.draw()

def update(frame):

    print(f'time: {frame}')
    ax.clear()
    data = probs_map[frame]
    # redraw faded fill + opaque lines

    artists['img'] = ax.contourf(
        lons, lats, probs_map[frame],
        levels=levels[1:], cmap=cmap, norm=norm,
        alpha=0.6, transform=ccrs.PlateCarree()
    )
    artists['cs'] = ax.contour(
        lons, lats, probs_map[frame],
        levels=levels[1:],
        colors=[cmap(norm(lv)) for lv in levels],
        linewidths=0.5,
        transform=ccrs.PlateCarree()
    )
    if frame>26:
        ds = Dataset(ens_files[0,frame], mode='r')
        refl_comp = get_reflectivity(ds)
        ds.close()

        refl_contour = ax.contour(lons, lats, refl_comp, levels=[15],
                colors='black', linewidths=1.5, transform=ccrs.PlateCarree())

    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='gray',linewidth=1)
    ax.coastlines()
    ax.set_extent([-104.5, -102.5, 40.5, 42])
    ax.add_feature(roads_feature, linewidth=0.45, edgecolor='darkgray')

    time_str = timestamps[frame].strftime('%Y-%m-%d %H%M')
    curr_str = f'Valid: {time_str} UTC'
    valid_time.set_text(curr_str)

    return artists['img'], artists['cs']

ani = animation.FuncAnimation(
        fig, update, frames=n_times, interval=200, blit=False
        )

ani.save('animation.gif', writer=PillowWriter(fps=5))
