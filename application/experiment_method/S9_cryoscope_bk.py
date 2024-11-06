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

from exp.save_data import save_nc, save_fig, create_folder
from exp.plotting import plot_and_save_cryoscope_bk
import matplotlib.pyplot as plt

# Set parameters
ro_elements = ["q0_ro", "q1_ro", "q2_ro", "q3_ro", "q4_ro"]
q_name = "q2_xy"
z_name = "q2_z"

n_avg = 2000
const_flux_len = 240 # unit ns <200
const_flux_amp =  0.22 #0.22, 0.18 ,0.145
pad_zeros = (20,0)
save_data = True
save_dir = link_config["path"]["output_root"]
save_name = f"{q_name}_cryoscope_bk"
folder_label = "cryoscope_bk_1" #your data and plots with be saved under a new folder with this name

from exp.cryoscope_bk import cryoscope_bk
dataset = cryoscope_bk( ro_elements, q_name, z_name, const_flux_amp, const_flux_len, n_avg, config, qmm, initializer=init_macro, pad_zeros=pad_zeros)

# "mixer", "r90", "time"
folder_save_dir = 0
if save_data: 
    folder_save_dir = create_folder(save_dir, folder_label)
    save_nc(folder_save_dir, save_name, dataset)

# Plot
plot_and_save_cryoscope_bk(dataset, folder_save_dir, pad_zeros, const_flux_len, save_data)