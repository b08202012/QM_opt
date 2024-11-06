import numpy as np
import matplotlib.pyplot as plt
from exp.cz_chavron import plot_cz_chavron, plot_cz_couplerz, plot_cz_ramsey
import xarray as xr
import os
filename = 'C:\\QM\\Data\\20240620_DR4_5Q4C_0430#7_new\\20241025_1752_q3q4_cz_chavron_shot.nc'
dataset = xr.open_dataset(filename,engine='netcdf4', format='NETCDF4')

time = dataset.coords["time"].values
#c_amps = dataset["c_amps"].values
amps = dataset.coords["amps"].values



threhold_q3 = 1.018e-4#2.615e-5#6.097e-5
threhold_q4 = 6.146e-5
population = np.zeros((len(amps),len(time)))
for i in range(500):
    for j in range(len(amps)):
        for k in range(len(time)):
            if dataset["q3_ro"].values[0,i,j,k] >= threhold_q3 and dataset["q4_ro"].values[0,i,j,k] >= threhold_q4:
                population[j,k] += 1
population = population/500
#print(population[10,10])
#print(amps)

#print(np.arange(-0.04,-0.05,-0.0001))
#fig, ax = plt.subplots()
#ax.plot(time,dataset["q1_ro"][0,11])

"""
z_offset = 0.05
z = 0.22
period = 0.65
Ec = 200
fc = 6.13 + (-221+0.75+5-1.25-0.5-10.89-0.264+1.255-0.48+1.03-0.06-0.23-1.1-0.8+0.63-1.76)/1000
phi = np.pi/period*(z-z_offset)
J_max = (fc*1000 +Ec)**2/(8*Ec*np.cos(phi))
zc = np.arange(0.2-0.17, 0.2+0.05, 0.001)
fq = (8*Ec*J_max*np.cos(np.pi/period*(c_amps-z_offset)))**0.5 - Ec
z_offset2 = 0.07
z2 = 0.148
period2 = 0.78
Ec2 = 200
fc2 = 4.7648+1.205e-3-0.073e-3
phi2 = np.pi/period2*(z2-z_offset2)
J_max2 = (fc2*1000 +Ec)**2/(8*Ec*(np.cos(phi2)**2+0.3**2*np.sin(phi2)**2)**0.5)
zc2 = np.arange(0.2-0.17, 0.2+0.05, 0.001)
fq2 = (8*Ec*J_max2*(np.cos(np.pi/period2*(amps-z_offset2))**2+0.3**2*np.sin(np.pi/period2*(amps-z_offset2))**2)**0.5)**0.5 - Ec
"""

for ro_name, data in dataset.data_vars.items():
    print(data.shape)
    fig, ax = plt.subplots()
    plot_cz_chavron(time,amps,population,ax)
    #plot_cz_couplerz(amps,c_amps,data.values[0],ax)
    #plot_cz_ramsey(c_amps, time, data.values[0],ax)
    #ax[1].plot(time, fq)
fig,ax = plt.subplots()
#ax.plot(time, dataset["q3_ro"].values[0,24])
#ax.set_xlabel("interaction time (ns)")
plt.show()
