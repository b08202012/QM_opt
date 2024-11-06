# Import necessary file
from pathlib import Path
link_path = Path(__file__).resolve().parent.parent/"config_api"/"config_link.toml"

from QM_driver_AS.ultitly.config_io import import_config, import_link
link_config = import_link(link_path)
config_obj, spec = import_config( link_path )

config = config_obj.get_config()
qmm, _ = spec.buildup_qmm()

from ab.QM_config_dynamic import initializer
#init_macro = initializer(200000,mode='wait')

from exp.save_data import save_nc, save_fig

import matplotlib.pyplot as plt
import numpy as np
# Set parameters
init_macro = initializer(100000,mode='wait')

ro_element=["q3_ro", "q4_ro"]
flux_Qi = 4
control_Qi = 4
target_Qi = 3
flux_Ci = 8 # coupler
preprocess = "shot"

n_avg = 1000
cz_amps_range = (-0.0325*2,-0.0275*2) #
cz_amps_resolution = 0.0001*2      #
cz_amps = 0.2105#0.0896#0.2105
cz_time_max = 100
cz_time_resolution = 4
cz_time = 152 # ns #
couplerz_amps_range = (-0.022*2,-0.012*2) #
couplerz_amps_resolution = 0.0002*2     #
c_amps = 0.0

save_data = True
save_dir = link_config["path"]["output_root"]
save_name = f"q{control_Qi}q{target_Qi}_cz_phasediff_1D"

from exp.cz_phase_diff import CZ_phase_diff, CZ_phase_compensate,CZ_phase_diff_time,CZ_phase_diff_1D
#dataset = CZ_phase_diff_1D(cz_amps,cz_time,ro_element,flux_Qi,control_Qi,target_Qi,flux_Ci,preprocess,qmm,config,n_avg=n_avg,initializer=init_macro,simulate=False)
#dataset = CZ_phase_diff(cz_amps_range,cz_amps_resolution,cz_time,couplerz_amps_range,couplerz_amps_resolution,ro_element,flux_Qi,control_Qi,target_Qi,flux_Ci,preprocess,qmm,config,n_avg=n_avg,initializer=init_macro,simulate=False)
#dataset = CZ_phase_diff_time(cz_amps_range,cz_amps_resolution,cz_time_max,cz_time_resolution,c_amps,ro_element,flux_Qi,control_Qi,target_Qi,flux_Ci,qmm,config,n_avg=n_avg,initializer=init_macro,simulate=False)
dataset = CZ_phase_compensate(c_amps,cz_amps,cz_time,ro_element,flux_Qi,target_Qi,flux_Ci,preprocess,qmm,config,n_avg=n_avg,initializer=init_macro,simulate=False)
if save_data: save_nc(save_dir, save_name, dataset)


