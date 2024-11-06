import xarray as xr
from exp.tomography import StateTomography, calculate_2Q_block_vector, calculate_block_vector, calculate_density_matrix, plot_density_matrix
import matplotlib.pyplot as plt


filename = 'C:\\QM\\Data\\20240620_DR4_5Q4C_0430#7_new\\20241024_1854_q3_roq4_ro_state_tomography.nc'
dataset = xr.open_dataset(filename,engine='netcdf4', format='NETCDF4')

threshold1 = 1.035e-4#6.840e-5#4.593e-5
threshold2 = 9.303e-5#1.434e-5#1.376e-4

data_i1 = dataset["q3_ro"].values[0].transpose((1,2,0))
data_i2 = dataset["q4_ro"].values[0].transpose((1,2,0))
print(data_i1.shape)
vect_dis = calculate_2Q_block_vector(data_i1, threshold1, data_i2, threshold2)
density_matrix = calculate_density_matrix(data_i1, threshold1, data_i2, threshold2)
print(density_matrix[0,0],density_matrix[1,1],density_matrix[2,2], density_matrix[3,3])
fig,ax = plt.subplots()
plot_density_matrix(density_matrix)
plt.show()