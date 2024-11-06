# Import necessary file
from pathlib import Path
link_path = Path(__file__).resolve().parent.parent/"config_api"/"config_link.toml"

from QM_driver_AS.ultitly.config_io import import_config, import_link
link_config = import_link(link_path)
config_obj, spec = import_config( link_path )

config = config_obj.get_config()
qmm, _ = spec.buildup_qmm()

from ab.QM_config_dynamic import initializer
init_macro = initializer(200000,mode='wait')

from exp.save_data import DataPackager

import matplotlib.pyplot as plt
import numpy as np
# Set parameters
init_macro = initializer(200000,mode='wait')

ro_element = ["q3_ro", "q4_ro"]
flux_Qi = 4
excited_Qi = [3,4]
flux_Ci = 8
n_avg = 500
preprocess = "shot" # ave or shot

time_max = 0.8 # us
time_resolution = 0.008 # us
z_amps = -0.03*2
z_amps_range = (0.1*2,0.12*2)
z_amps_resolution = 0.0002*2
coupler_z = 0.0*2
couplerz_amps_range = (0.0*2,0.3*2)
couplerz_amps_resolution = 0.003*2

save_data = True
save_dir = link_config["path"]["output_root"]
save_name = f"q{excited_Qi[0]}q{excited_Qi[1]}_cz_chavron_shot"

from exp.cz_chavron import CZ,CZ_couplerz, CZ_coupler_time
# dataset = CZ(time_max,time_resolution,z_amps_range,z_amps_resolution,ro_element,flux_Qi,excited_Qi,flux_Ci,coupler_z,preprocess,qmm,config,n_avg=n_avg,initializer=init_macro,simulate=False)
dataset = CZ_coupler_time(time_max,time_resolution,couplerz_amps_range,couplerz_amps_resolution,ro_element,flux_Qi,excited_Qi,flux_Ci,qubit_z,preprocess,qmm,config,n_avg=n_avg,initializer=init_macro,simulate=False)
# dataset = CZ_couplerz(z_amps_range,z_amps_resolution,couplerz_amps_range,couplerz_amps_resolution,ro_element,flux_Qi,excited_Qi,flux_Ci,preprocess,qmm,config,n_avg=n_avg,initializer=init_macro,simulate=False)
folder_label = "CZ" #your data and plots will be saved under a new folder with this name
save_data = 1
if save_data: 
    from exp.save_data import DataPackager
    save_dir = link_config["path"]["output_root"]
    dp = DataPackager( save_dir, folder_label )
    dp.save_config(config)
    dp.save_nc(dataset,save_name) 

# Plot
save_figure = 1
# import xarray as xr
# import matplotlib.pyplot as plt
# dataset = xr.open_dataset(r"C:\Users\quant\SynologyDrive\09 Data\Fridge Data\Qubit\20240920_DRKe_5Q4C\save_data\CZ_sweet\crosstalk_not_compensated\20241004_055553_CZ\q1q0_cz_couplerz.nc")

time = dataset.coords["time"].values
# coupler_flux = dataset.coords["c_amps"].values
flux = dataset.coords["amps"].values

import numpy as np
from exp.cz_chavron import plot_cz_chavron,plot_cz_couplerz, plot_coupler_z_vs_time, plot_cz_frequency_vs_flux, plot_cz_period_vs_flux
for ro_name, data in dataset.data_vars.items():
    fig, ax = plt.subplots()
    # plot_cz_chavron(time,flux,data.values[0],ax)
    # plot_cz_frequency_vs_flux(time, flux, data.values[0], np.inf,ax)
    # plot_cz_period_vs_flux(time, flux, data.values[0], np.inf,ax)
    plot_coupler_z_vs_time(time, flux, data.values[0],ax)
    # plot_cz_couplerz(flux,coupler_flux,data.values[0],ax)
    # if save_figure: dp.save_fig(fig, f"{save_name}_{ro_name}")
plt.show()
