import numpy as np
import matplotlib.pyplot as plt
from exp.cz_chavron import plot_cz_chavron
import xarray as xr
import os
filename = 'C:\\QM\\Data\\20240620_DR4_5Q4C_0430#7\\20240813_1701_q4_xy_piscope_cc.nc'
dataset = xr.open_dataset(filename,engine='netcdf4', format='NETCDF4')

qubit_freq = 5.0505+0.281e-3-0.101e-3+0.365e-3

time = dataset.coords["pi_timing"].values
freq = qubit_freq*1000+dataset.coords["frequency"].values

for ro_name, data in dataset.data_vars.items():
    fig_0, ax_0 = plt.subplots()
    fig, ax = plt.subplots()
    ax.set_title('init=100us zpulse_time=1000+10us')
    ax.set_xlabel("driving freq (MHz)")
    ax.set_ylabel("time (ns)")
    pcm = ax.pcolormesh( freq, time, data.values[0], cmap='RdBu')# , vmin=z_min, vmax=z_max)
    plt.colorbar(pcm, label='Value')
plt.show()