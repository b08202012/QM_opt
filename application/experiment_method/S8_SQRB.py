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

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

# Set parameters
init_macro = initializer(300000,mode='wait')
# Set parameters
ro_elements = ["q3_ro"]
q_name = "q3_xy"

gate_length = 40

num_of_sequences = 100  # Number of random sequences
n_avg = 20  # Number of averaging loops for each random sequence
max_circuit_depth = 200  # Maximum circuit depth
delta_clifford = 5  #  Play each sequence with a depth step equals to 'delta_clifford - Must be > 1
assert (max_circuit_depth / delta_clifford).is_integer(), "max_circuit_depth / delta_clifford must be an integer."
seed = 345324  # Pseudo-random number generator seed

save_data = True
save_dir = link_config["path"]["output_root"]
save_name = f"{q_name[0]}_sqrb"
folder_label = "SQRB_1" #your data and plots with be saved under a new folder with this name

from exp.randomized_banchmarking_sq import single_qubit_RB
x, value_avg, error_avg = single_qubit_RB( num_of_sequences, max_circuit_depth, delta_clifford, q_name, ro_elements, config, qmm, seed=seed, gate_length=gate_length)
print(x, value_avg)
output_data = {}
output_data[q_name] = (["gate_number"], np.array(value_avg)) 
output_data[f"{q_name}_err"] = (["gate_number"], np.array(error_avg))

dataset = xr.Dataset(
    output_data,
    coords={"gate_number": x, }
)    

folder_save_dir = 0
if save_data: 
    folder_save_dir = create_folder(save_dir, folder_label)
    save_nc(folder_save_dir, save_name, dataset)

# Plot
from exp.randomized_banchmarking_sq import plot_SQRB_result
fig, ax = plt.subplots(2)

plot_SQRB_result( x, value_avg, error_avg )
if save_data: save_fig(folder_save_dir, save_name, dataset)

plt.show()