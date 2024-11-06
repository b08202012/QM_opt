import numpy as np
import matplotlib.pyplot as plt

# Dynamic config
# Import necessary file
from pathlib import Path
link_path = Path(__file__).resolve().parent.parent/"config_api"/"config_link.toml"

from QM_driver_AS.ultitly.config_io import import_config, import_link
link_config = import_link(link_path)
config_obj, spec = import_config( link_path )

config = config_obj.get_config()
qmm, _ = spec.buildup_qmm()

# Dynamic config
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
link_config = import_link(link_path)


ro_elements = ['q0_ro','q1_ro']
excited_q = 'q1_xy'
z_name = ['q1_z']


n_avg = 200  # The number of averages

# amps = np.arange(0.10, 0.14, 0.002)  # The relative flux amplitude absZvolt-offset
# time = np.arange(40, 8000, 200) # The flux pulse durations in clock cycles (4ns) - Must be larger than 4 clock cycles.

mid = 0.18
ra = 0.05
amps = np.arange(mid-ra, mid+ra, ra/25)  # The relative flux amplitude absZvolt-offset
time = np.arange(16, 1000, 4) # The flux pulse durations in clock cycles (4ns) - Must be larger than 4 clock cycles.
coupler_z = "q2_z"
coupler_amp = +0.0
cc = time/4

from exp.iSWAP_J import exp_coarse_iSWAP, plot_ana_iSWAP_chavron
dataset = exp_coarse_iSWAP( coupler_z, coupler_amp, amps, cc, excited_q, ro_elements, z_name, config, qmm, n_avg=n_avg, simulate=False, initializer=init_macro )

save_data = True
save_dir = link_config["path"]["output_root"]
save_name = f"iSWAP"
folder_label = "find_iSWAP_1"
if save_data: 
    save_dir = create_folder(save_dir, folder_label)
    save_nc(save_dir, save_name, dataset)

time = dataset.coords["time"].values*4
amps = dataset.coords["amplitude"].values
for ro_name, data in dataset.data_vars.items():

    fig, ax = plt.subplots(2)
    plot_ana_iSWAP_chavron( data.values, amps, time, ax )
    ax[0].set_title(ro_name)
    ax[1].set_title(ro_name)

plt.show()



   