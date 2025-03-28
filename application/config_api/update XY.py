

from pathlib import Path
link_path = Path(__file__).resolve().parent/"config_link.toml"

from QM_driver_AS.ultitly.config_io import import_config
config_obj, spec = import_config( link_path )

from qspec.update import update_controlFreq, update_controlWaveform

import numpy as np

xy_infos = {

# e.g.
# "q0":{   
#     "q_freq": 3.05,   # GHz
#     "LO": 3.0,        # GHz
#     "pi_amp": 0.15,   # factor
#     "pi_len": 40,     # ns 
#     "90_corr": 1
# },

    "q0":{   # unfinish
        "q_freq": 3.05,     # GHz
        "LO": 3.0,          # GHz
        "pi_amp": 0.15,
        "pi_len": 40,
        "90_corr": 1
    },
    "q1":{   
        "q_freq": 4.049,    # GHz
        "LO": 4.2,          # GHz
        "pi_amp": 0.15*0.58,
        "pi_len": 40,
        "90_corr": 1
    },
    "q2":{   # unfinish, maybe low T1
        "q_freq": 4.79, # GHz
        "LO": 4.7, # GHz
        "pi_amp": 0.5,
        "pi_len": 1000,
        "90_corr": 1,
        # "wf": "sin",
    },
    "q3":{   
        "q_freq": 3.627, # GHz
        "LO": 3.7, # GHz
        "pi_amp": 0.5,
        "pi_len": 40,
        "90_corr": 1
    },
    "q3":{
        "q_freq": 4.69966,#4.6665+2.349e-3-1.621e-3+0.332e-3-1.126e-3-0.764e-3-2.067e-3-0.547e-3+0.074e-3,#4.6689-2.897e-3+0.415e-3+2.402e-3-1.236e-3,#4.546, # GHz
        "LO": 4.85, # GHz
        "pi_amp": 0.1114941471,
        "pi_len": 40,
        "90_corr": 1.01578436269,
        'drag_coef': 0.3,
        "anharmonicity": -200,
        "AC_stark_detuning" : 0,
        "wf": {
            'dragg':{
                'sfactor' : 4,
            }
        }
    },
    "q5":{   # unfinish, maybe low T1
        "q_freq": 4.45, # GHz
        "LO": 4.4, # GHz
        "pi_amp": 0.5,
        "pi_len": 40,
        "90_corr": 1
    },
}    
updating_qubit = ["q0","q3"]

# name, q_freq(GHz), LO(GHz), amp, len, half
# update_info = [['q1', 4.526-0.005, 4.60, 0.2*0.9*1.05*1.01, 40, 1.03 ]]
for i in updating_qubit:
    # wiring = spec.get_spec_forConfig('wire')
    q_name = i
    qubit_LO = xy_infos[i]["LO"]
    qubit_RF = xy_infos[i]["q_freq"]
    ref_IF = (qubit_RF-qubit_LO)*1000
    print(f"center {ref_IF} MHz")

    pi_amp = xy_infos[i]["pi_amp"] if "pi_amp" in xy_infos[i].keys() else 0.1
    pi_len = xy_infos[i]["pi_len"] if "pi_len" in xy_infos[i].keys() else 40
    corr_90 = xy_infos[i]["90_corr"] if "90_corr" in xy_infos[i].keys() else 1
    drag_co = xy_infos[i]["drag_coef"] if "drag_coef" in xy_infos[i].keys() else 0.5
    anharm = xy_infos[i]["anharmonicity"] if "anharmonicity" in xy_infos[i].keys() else -200
    AC_stark = xy_infos[i]["AC_stark_detuning"] if "AC_stark_detuning" in xy_infos[i].keys() else 0
    qubit_wf = xy_infos[i]["wf"] if "wf" in xy_infos[i].keys() else 0

    updated_xy_freq = spec.update_aXyInfo_for(target_q=q_name, 
                                              IF=ref_IF,
                                              LO=qubit_LO,
                                              amp=pi_amp, 
                                              len=pi_len, 
                                              half=corr_90,
                                              draga=drag_co,
                                              anh=anharm,
                                              ac=AC_stark,
                                              wf=qubit_wf
                                              )
    
    update_controlFreq(config_obj, updated_xy_freq)
    if np.abs(ref_IF) > 350:
        print("Warning IF > +/-350 MHz, IF is set 350 MHz")
        ref_IF = np.sign(ref_IF)*350

    update_controlWaveform(config_obj, spec.get_spec_forConfig("xy"), target_q=q_name )

from QM_driver_AS.ultitly.config_io import output_config
output_config( link_path, config_obj, spec )
