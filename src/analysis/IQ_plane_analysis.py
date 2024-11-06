import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

filename = 'C:\\QM\\Data\\20240620_DR4_5Q4C_0430#7_new\\20241025_1513_q3_roq4_ro_state_tomography.nc'
dataset = xr.open_dataset(filename,engine='netcdf4', format='NETCDF4')

threhold_q3 = 1.083e-4#6.571e-5#4.77e-5#2.615e-5#6.097e-5
threhold_q4 = 6.146e-5

q3_excite_I = []
q3_excite_Q = []
q3_ground_I = []
q3_ground_Q = []
q4_excite_I = []
q4_excite_Q = []
q4_ground_I = []
q4_ground_Q = []
q3_0_q4_2_I = []
q3_0_q4_2_Q = []
shot_num = 1000
fig,ax = plt.subplots(2,2)
for i in range(shot_num):
    if dataset["q3_ro"].values[0,i,0,0] >= threhold_q3:
        q3_excite_I.append(dataset["q3_ro"].values[0,i,0,0])
        q3_excite_Q.append(dataset["q3_ro"].values[1,i,0,0])
    else:
        q3_ground_I.append(dataset["q3_ro"].values[0,i,0,0])
        q3_ground_Q.append(dataset["q3_ro"].values[1,i,0,0])
    if dataset["q4_ro"].values[0,i,0,0] >= threhold_q4:
        q4_excite_I.append(dataset["q4_ro"].values[0,i,0,0])
        q4_excite_Q.append(dataset["q4_ro"].values[1,i,0,0])
    else:
        q4_ground_I.append(dataset["q4_ro"].values[0,i,0,0])
        q4_ground_Q.append(dataset["q4_ro"].values[1,i,0,0])
    if dataset["q3_ro"].values[0,i,0,0] < threhold_q3 and dataset["q4_ro"].values[0,i,0,0] >= threhold_q4:
        q3_0_q4_2_I.append(dataset["q4_ro"].values[0,i,0,0])
        q3_0_q4_2_Q.append(dataset["q4_ro"].values[1,i,0,0])
ax[0,0].plot(q3_ground_I, q3_ground_Q, ".", label="Ground", markersize=2)
ax[0,0].plot(q3_excite_I, q3_excite_Q, ".", label="Excite", markersize=2)
ax[0,0].axvline(x=threhold_q3, color="k", ls="--", alpha=0.5)
ax[0,1].plot(dataset["q3_ro"].values[0,:,0,0], dataset["q3_ro"].values[1,:,0,0], ".", markersize=2)
ax[1,0].plot(q4_ground_I, q4_ground_Q, ".", label="Ground", markersize=2)
ax[1,0].plot(q4_excite_I, q4_excite_Q, ".", label="Excite", markersize=2)
ax[1,0].plot(q3_0_q4_2_I, q3_0_q4_2_Q, ".", label="Second excite", markersize=2)
ax[1,1].plot(dataset["q4_ro"].values[0,:,0,0], dataset["q4_ro"].values[1,:,0,0], ".", markersize=2)
ax[1,0].axvline(x=threhold_q4, color="k", ls="--", alpha=0.5)
plt.legend()
plt.show()