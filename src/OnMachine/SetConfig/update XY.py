from OnMachine.SetConfig.config_path import spec_loca, config_loca
from config_component.configuration import import_config, configuration_read_dict
from config_component.channel_info import import_spec


spec = import_spec( spec_loca )
config_obj = import_config( config_loca )

from config_component.update import update_controlFreq, update_controlWaveform

import numpy as np

xy_infos = [
    {
        "name":"q1",
        "q_freq": 4.550, # GHz
        "LO": 4.6, # GHz
        "pi_amp": 0.2*0.9*1.05*1.01*0.7*0.75,
        "pi_len": 40,
        "90_corr": 1.03
    },
    {
        "name":"q2",
        "q_freq": 4.2, # GHz
        "LO": 4.11, # GHz
        "pi_amp": 0.3*1.1*1.02*0.7*1.2*1.02,
        "pi_len": 40,
        "90_corr": 1.2*0.92
    },
    {
        "name":"q5",
        "q_freq": 5.0, # GHz
        "LO": 5.0, # GHz
        "pi_amp": 0.2*0.9*1.05*1.01*0.7,
        "pi_len": 40,
        "90_corr": 1.03
    }
]
# name, q_freq(GHz), LO(GHz), amp, len, half
# update_info = [['q1', 4.526-0.005, 4.60, 0.2*0.9*1.05*1.01, 40, 1.03 ]]
for i in xy_infos:
    # wiring = spec.get_spec_forConfig('wire')
    q_name = i["name"]
    qubit_LO = i["LO"]
    qubit_RF = i["q_freq"]
    ref_IF = (qubit_RF-qubit_LO)*1000

    print(f"center {ref_IF} MHz")
    pi_amp = i["pi_amp"]
    pi_len = i["pi_len"]

    update_controlFreq(config_obj, spec.update_aXyInfo_for(target_q=q_name,IF=ref_IF,LO=qubit_LO))
    if np.abs(ref_IF) > 350:
        print("Warning IF > +/-350 MHz, IF is set 350 MHz")
        ref_IF = np.sign(ref_IF)*350

    spec.update_aXyInfo_for(target_q=q_name, amp=pi_amp, len=pi_len, half=i["90_corr"])
    update_controlWaveform(config_obj, spec.get_spec_forConfig("xy"), target_q=q_name )

import json
file_path = 'output.json'
# Open the file in write mode and use json.dump() to export the dictionary to JSON
with open(file_path, 'w') as json_file:
    json.dump(config_obj.get_config(), json_file, indent=2)

spec.export_spec(spec_loca)
config_obj.export_config(config_loca)