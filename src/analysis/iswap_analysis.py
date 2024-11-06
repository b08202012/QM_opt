import numpy as np
import matplotlib.pyplot as plt
from exp.cz_chavron import plot_cz_chavron, plot_cz_couplerz
import xarray as xr
from exp.iSWAP_J import exp_coarse_iSWAP, plot_ana_iSWAP_chavron
import os
filename = 'C:\\QM\\Data\\20240916_DR4_2Q2C_0719\\20240919_1503_iswap_q1_xy_q0_z_q2_z_0.0.nc'
dataset = xr.open_dataset(filename,engine='netcdf4', format='NETCDF4')

time = dataset.coords["time"].values*4
amps = dataset.coords["amplitude"].values
for ro_name, data in dataset.data_vars.items():

    fig, ax = plt.subplots(2)
    plot_ana_iSWAP_chavron( data.values, amps, time, ax )
    ax[0].set_title(ro_name)
    ax[1].set_title(ro_name)
    fig, ax = plt.subplots()
    ax.plot(time,data.values[0,:,11])
plt.show()