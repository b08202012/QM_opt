import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

filename = 'C:\\QM\\Data\\20240620_DR4_5Q4C_0430#7_new\\20241025_1519_q3q4_cz_phasediff_1D.nc'
dataset = xr.open_dataset(filename,engine='netcdf4', format='NETCDF4')
cz_amp = [1]#dataset.coords["cz_amp"].values
c_amp = [1]#dataset.coords["c_amp"].values
#cz_time = dataset.coords["cz_time"].values

threhold_q3 = 1.083e-4#6.571e-5#4.77e-5#2.615e-5#6.097e-5
threhold_q4 = 6.146e-5#1.079e-5#1.303e-4#-1.483e-6#8.246e-7
ro_element = "q4_ro" # target qubit readout
shot_num = 1000
data = dataset[ro_element].values[0]
print(data.shape)

population_gy = np.zeros((len(c_amp),len(cz_amp))) # ground state y-basis
population_gx = np.zeros((len(c_amp),len(cz_amp))) # ground state x-basis
population_gz = np.zeros((len(c_amp),len(cz_amp))) # ground state z-basis

population_ey = np.zeros((len(c_amp),len(cz_amp))) # excited state y-basis
population_ex = np.zeros((len(c_amp),len(cz_amp))) # excited state x-basis
population_ez = np.zeros((len(c_amp),len(cz_amp))) # excited state z-basis


phase_map = np.zeros((len(c_amp),len(cz_amp)))
phase_map_small = np.zeros((len(c_amp),len(cz_amp)))

e_sqrt_map = np.zeros((len(c_amp),len(cz_amp)))
e_Sz_map = np.zeros((len(c_amp),len(cz_amp)))
g_sqrt_map = np.zeros((len(c_amp),len(cz_amp)))
g_Sz_map = np.zeros((len(c_amp),len(cz_amp)))



for i in range(shot_num):
    for j in range(len(c_amp)):
        for k in range(len(cz_amp)):

            if data[i,j,k,0,0] >= threhold_q3:
                population_gx[j,k] += 1
            if data[i,j,k,0,1] >= threhold_q3:
                population_gy[j,k] += 1
            if data[i,j,k,0,2] >= threhold_q3:
                population_gz[j,k] += 1
            if data[i,j,k,1,0] >= threhold_q3:
                population_ex[j,k] += 1
            if data[i,j,k,1,1] >= threhold_q3:
                population_ey[j,k] += 1
            if data[i,j,k,1,2] >= threhold_q3:
                population_ez[j,k] += 1
# a=87.3/(87.3+11.5)
population_gy = population_gy/shot_num#*a
population_gx = population_gx/shot_num#*a
population_gz = population_gz/shot_num#*a

population_ey = population_ey/shot_num#*a
population_ex = population_ex/shot_num#*a
population_ez = population_ez/shot_num#*a


for i in range(len(c_amp)):
    for j in range(len(cz_amp)):
        gy = (1-population_gy[i,j])-population_gy[i,j]
        gx = (1-population_gx[i,j])-population_gx[i,j]
        gz = (1-population_gz[i,j])-population_gz[i,j] 
        g_phi = np.arctan2(gy,gx)

        ey = (1-population_ey[i,j])-population_ey[i,j]
        ex = (1-population_ex[i,j])-population_ex[i,j]
        ez = (1-population_ez[i,j])-population_ez[i,j] 
        e_phi = np.arctan2(ey,ex)

        if abs((abs(g_phi-e_phi)-np.pi))<=0.1*2*np.pi:
            phase_map[i,j] = (abs(g_phi-e_phi)-np.pi)
            phase_map_small[i,j] = (abs(g_phi-e_phi)-np.pi)

        else:
            phase_map[i,j] = (abs(g_phi-e_phi)-np.pi)
            phase_map_small[i,j] = 0.1*2*np.pi
        # if abs(abs(g_phi-e_phi)-np.pi) <0.01:
        #     print(c_amp[i],cz_amp[j])
        #     print((g_phi-e_phi))
        #     print(g_phi,e_phi)
        #     print(gy,gx)
        #     print(gy**2+gx**2)
        #     print(ey,ex)
        #     print(ey**2+ex**2)
        #     print(population_gy[i,j])
        #     print("------------")
        g_sqrt_map[i, j] = np.sqrt(gx**2 + gy**2)
        g_Sz_map[i, j] = gz
        e_sqrt_map[i, j] = np.sqrt(ex**2 + ey**2)
        e_Sz_map[i, j] = ez
        #print((g_phi-e_phi))


"""
z_offset = 0.05
z = 0.22
period = 0.65
Ec = 200
fc = 6.13 + (-221+0.75+5-1.25-0.5-10.89-0.264+1.255-0.48+1.03-0.06-0.23-1.1-0.8+0.63-1.76)/1000
phi = np.pi/period*(z-z_offset)
J_max = (fc*1000 +Ec)**2/(8*Ec*np.cos(phi))
zc = np.arange(0.2-0.17, 0.2+0.05, 0.001)
fq = (8*Ec*J_max*np.cos(np.pi/period*(c_amp-z_offset)))**0.5 - Ec
z_offset2 = 0.07
z2 = 0.148
period2 = 0.78
Ec2 = 200
fc2 = 4.7648+1.205e-3-0.073e-3
phi2 = np.pi/period2*(z2-z_offset2)
J_max2 = (fc2*1000 +Ec)**2/(8*Ec*(np.cos(phi2)**2+0.3**2*np.sin(phi2)**2)**0.5)
zc2 = np.arange(0.2-0.17, 0.2+0.05, 0.001)
fq2 = (8*Ec*J_max2*(np.cos(np.pi/period2*(cz_amp-z_offset2))**2+0.3**2*np.sin(np.pi/period2*(cz_amp-z_offset2))**2)**0.5)**0.5 - Ec
"""


# 第一張圖
fig1, ax1 = plt.subplots()
a1 = ax1.pcolormesh(cz_amp, c_amp, phase_map, cmap='RdBu')
ax1.set_xlabel("q2 z")
ax1.set_ylabel("coupler z")
plt.colorbar(a1, ax=ax1, label="phase difference - pi")

# 第二張圖
fig2, ax2 = plt.subplots()
a2 = ax2.pcolormesh(cz_amp, c_amp, e_sqrt_map, cmap='RdBu')  # 替換phase_map或其他資料
ax2.set_xlabel("q2 z")
ax2.set_ylabel("coupler z")
plt.colorbar(a2, ax=ax2, label="e_sqrt(Sx\^2+Sy\^2)")

# 第三張圖
fig3, ax3 = plt.subplots()
a3 = ax3.pcolormesh(cz_amp, c_amp, e_Sz_map, cmap='RdBu')  # 替換phase_map或其他資料
ax3.set_xlabel("q2 z")
ax3.set_ylabel("coupler z")
plt.colorbar(a3, ax=ax3, label="e_Sz")

# 第四張圖
fig4, ax4 = plt.subplots()
a4 = ax4.pcolormesh(cz_amp, c_amp, g_sqrt_map, cmap='RdBu')  # 替換phase_map或其他資料
ax4.set_xlabel("q2 z")
ax4.set_ylabel("coupler z")
plt.colorbar(a4, ax=ax4, label="g_sqrt(Sx\^2+Sy\^2)")

# 第五張圖
fig5, ax5 = plt.subplots()
a5 = ax5.pcolormesh(cz_amp, c_amp, g_Sz_map, cmap='RdBu')  # 替換phase_map或其他資料
ax5.set_xlabel("q2 z")
ax5.set_ylabel("coupler z")
plt.colorbar(a5, ax=ax5, label="g_Sz")

# 第一張圖
fig1, ax6 = plt.subplots()
a6 = ax6.pcolormesh(cz_amp, c_amp, phase_map_small, cmap='RdBu')
ax6.set_xlabel("q2 z")
ax6.set_ylabel("coupler z")
plt.colorbar(a6, ax=ax6, label="phase difference - pi")

plt.show()

