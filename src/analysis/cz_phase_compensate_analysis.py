import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

filename = 'C:\\QM\\Data\\20240620_DR4_5Q4C_0430#7_new\\20241025_1523_q4q3_cz_phasediff_1D.nc'
dataset = xr.open_dataset(filename,engine='netcdf4', format='NETCDF4')


threhold = 1.018e-4#1.035e-4#9.303e-5#6.571e-5#1.079e-5#8.246e-7
ro_element = "q3_ro" # target qubit readout
shot_num = 1000
data = dataset[ro_element].values[0]
print(data.shape)

population_y = 0
population_x = 0
for i in range(shot_num):
    if data[i,0] < threhold:
        population_y += 1
    if data[i,1] < threhold:
        population_x += 1
population_y = population_y/shot_num
population_x = population_x/shot_num

y = 2*population_y - 1
x = 2*population_x - 1
phi = np.arctan2(y,x)

print(phi)