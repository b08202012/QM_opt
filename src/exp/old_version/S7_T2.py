# Import necessary file
from pathlib import Path
link_path = Path(__file__).resolve().parent.parent/"config_api"/"config_link.toml"

from QM_driver_AS.ultitly.config_io import import_config, import_link
link_config = import_link(link_path)
config_obj, spec = import_config( link_path )

config = config_obj.get_config()
qmm, _ = spec.buildup_qmm()

from ab.QM_config_dynamic import initializer

from exp.save_data import save_nc, save_fig

import matplotlib.pyplot as plt

# Set parameters
init_macro = initializer(100000,mode='wait')
ro_elements = ["q0_ro", "q1_ro", "q2_ro", "q3_ro", "q4_ro"]
q_name = ["q4_xy"]


n_avg = 500
virtual_detune = 5

time_max = 10.0
time_resolution = 0.02
save_data = True
save_dir = link_config["path"]["output_root"]
save_name = f"{q_name[0]}_T2"

from exp.ramsey import exp_ramsey, multi_T2_exp
dataset = exp_ramsey( time_max, time_resolution, ro_elements, q_name, n_avg, config, qmm, virtual_detune=virtual_detune, initializer=init_macro)

repeat = 10
#dataset = multi_T2_exp( repeat, time_max, time_resolution, ro_elements, q_name, n_avg, config, qmm, virtual_detune=virtual_detune, initializer=init_macro)

if save_data: save_nc(save_dir, save_name, dataset)

# Plot
from exp.ramsey import plot_ramsey_oscillation, plot_multiT2
time = dataset.coords["time"].values
for ro_name, data in dataset.data_vars.items():
    fig, ax = plt.subplots(2)

    plot_ramsey_oscillation(time, data[0], ax[0])
    plot_ramsey_oscillation(time, data[1], ax[1])
    #rep = dataset.coords["repetition"].values
    #plot_multiT2( data, rep, time )
if save_data: save_fig(save_dir, save_name, dataset)

plt.show()

from exp.repetition_measurement import RepetitionMeasurement
re_exp = RepetitionMeasurement()
re_exp.exp_list = [my_exp]
re_exp.exp_name = ["T2"]
my_exp.shot_num = 400
dataset = re_exp.run(50)
save_name = f"{my_exp.xy_elements[0]}_EchoT2_stat"
if save_data: save_nc(save_dir, save_name, dataset["spin_echo"])