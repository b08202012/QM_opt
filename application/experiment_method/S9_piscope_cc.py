# Import necessary file
from pathlib import Path
link_path = Path(__file__).resolve().parent.parent/"config_api"/"config_link.toml"

from QM_driver_AS.ultitly.config_io import import_config, import_link
link_config = import_link(link_path)
config_obj, spec = import_config( link_path )

config = config_obj.get_config()
qmm, _ = spec.buildup_qmm()

from ab.QM_config_dynamic import initializer

from exp.save_data import save_nc, save_fig, create_folder
from exp.plotting import plot_and_save_piscope

# Set parameters

from exp.detuned_rabi_flux_pulse import DetunedRabiFluxPulse
my_exp = DetunedRabiFluxPulse(config, qmm)
my_exp.initializer = initializer(100000,mode='wait')
my_exp.ro_elements = ["q3_ro","q4_ro"]
my_exp.xy_elements = ["q4_xy"]
my_exp.z_elements = ["q4_z"]
my_exp.freq_range = (-50, 150)
my_exp.freq_resolution =2.0
my_exp.duration = 160
my_exp.pad_zeros = (80, 80) #ns
my_exp.time_resolution = 4 #ns
my_exp.amp_modify = -0.062 #x0.5 is voltage     0.44 0.36 0.29
 
dataset = my_exp.run(200)

save_data = True
save_dir = link_config["path"]["output_root"]
save_name = f"{my_exp.xy_elements[0]}_piscope_cc"
folder_label = "piscope_1" #your data and plots will be saved under a new folder with this name

if save_data: 
    save_dir = create_folder(save_dir, folder_label)
    save_nc(save_dir, save_name, dataset)

# Plot
plot_and_save_piscope(dataset, save_dir, save_data)
