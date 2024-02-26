"""
        CZ CHEVRON - 4ns granularity
The goal of this protocol is to find the parameters of the CZ gate between two flux-tunable qubits.
The protocol consists in flux tuning one qubit to bring the |11> state on resonance with |20>.
The two qubits must start in their excited states so that, when |11> and |20> are on resonance, the state |11> will
start acquiring a global phase when varying the flux pulse duration.

By scanning the flux pulse amplitude and duration, the CZ chevron can be obtained and post-processed to extract the
CZ gate parameters corresponding to a single oscillation period such that |11> pick up an overall phase of pi (flux
pulse amplitude and interation time).

This version sweeps the flux pulse duration using real-time QUA, which means that the flux pulse can be arbitrarily long
but the step must be larger than 1 clock cycle (4ns) and the minimum pulse duration is 4 clock cycles (16ns).

Prerequisites:
    - Having found the resonance frequency of the resonator coupled to the qubit under study (resonator_spectroscopy).
    - Having found the qubits maximum frequency point (qubit_spectroscopy_vs_flux).
    - Having calibrated qubit gates (x180) by running qubit spectroscopy, rabi_chevron, power_rabi, Ramsey and updated the configuration.
    - (Optional) having corrected the flux line distortions by running the Cryoscope protocol and updating the filter taps in the configuration.

Next steps before going to the next node:
    - Update the CZ gate parameters in the configuration.
"""

from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm import SimulationConfig
from configuration import *
import matplotlib.pyplot as plt
from qualang_tools.loops import from_array
from qualang_tools.results import fetching_tool
from qualang_tools.plot import interrupt_on_close
from qualang_tools.results import progress_counter
from common_fitting_func import *
import numpy as np
from macros import qua_declaration, multiplexed_readout
import warnings

warnings.filterwarnings("ignore")


####################
# Define variables #
####################
qubit_to_flux_tune = 4  # Qubit number to flux-tune
scale_reference = const_flux_amp # for const # 0.45
# scale_reference = cz_point_1_2_q2-idle_q2 # for gaussian-like cz

Qi = 4
n_avg = 100  # The number of averages
ts = np.arange(4, 60, 1)  # The flux pulse durations in clock cycles (4ns) - Must be larger than 4 clock cycles.
operation_flux_point = [0, -3.000e-01, -0.2525, -0.3433, -3.400e-01] 
amps = (np.arange(0.11, 0.56, 0.003)) 
q_id = [1,2,3,4]

res_F3 = cosine_func( operation_flux_point[2], *g1[2])
res_IF3 = (res_F3 - resonator_LO)/1e6
res_IF3 = int(res_IF3 * u.MHz)

res_F4 = cosine_func( operation_flux_point[3], *g1[3])
res_IF4 = (res_F4 - resonator_LO)/1e6
res_IF4 = int(res_IF4 * u.MHz)

with program() as cz:
    I, I_st, Q, Q_st, n, n_st = qua_declaration(nb_of_qubits=4)
    t = declare(int)  # QUA variable for the flux pulse duration
    a = declare(fixed)  # QUA variable for the flux pulse amplitude pre-factor.
    resonator_freq3 = declare(int, value=res_IF3)
    resonator_freq4 = declare(int, value=res_IF4)
    for i in q_id:
        set_dc_offset("q%s_z"%(i+1), "single", operation_flux_point[i])
    wait(flux_settle_time * u.ns)
    update_frequency(f"rr{3}", resonator_freq3)
    update_frequency(f"rr{4}", resonator_freq4)
    with for_(n, 0, n < n_avg, n + 1):
        with for_(*from_array(t, ts)):
            with for_(*from_array(a, amps)):
                # Put the two qubits in their excited states
                # play("x180", "q3_xy")
                play("x180", "q4_xy")
                align()
                # Wait some time to ensure that the flux pulse will arrive after the x90 pulse
                wait(flux_settle_time * u.ns)
                # Play a flux pulse on the qubit with the highest frequency to bring it close to the excited qubit while
                # varying its amplitude and duration in order to observe the SWAP chevron.
                play("const" * amp(a), f"q{qubit_to_flux_tune}_z", duration=t)
                # play("cz_1_2" * amp(a), f"q{qubit_to_flux_tune}_z", duration=t)
                align()
                wait(flux_settle_time * u.ns)
                align()
                multiplexed_readout(I, I_st, Q, Q_st, resonators=[1, 2, 3, 4], weights="rotated_")
                wait(thermalization_time * u.ns)
        save(n, n_st)

    with stream_processing():
        # for the progress counter
        n_st.save("n")
        # resonator 1
        I_st[0].buffer(len(amps)).buffer(len(ts)).average().save("I1")
        Q_st[0].buffer(len(amps)).buffer(len(ts)).average().save("Q1")
        # resonator 2
        I_st[1].buffer(len(amps)).buffer(len(ts)).average().save("I2")
        Q_st[1].buffer(len(amps)).buffer(len(ts)).average().save("Q2")
        # resonator 3
        I_st[2].buffer(len(amps)).buffer(len(ts)).average().save("I3")
        Q_st[2].buffer(len(amps)).buffer(len(ts)).average().save("Q3")
        # resonator 4
        I_st[3].buffer(len(amps)).buffer(len(ts)).average().save("I4")
        Q_st[3].buffer(len(amps)).buffer(len(ts)).average().save("Q4")



#####################################
#  Open Communication with the QOP  #
#####################################

qmm = QuantumMachinesManager(host=qop_ip, port=qop_port, cluster_name=cluster_name, octave=octave_config)

###########################
# Run or Simulate Program #
###########################
simulate = False

if simulate:
    # Simulates the QUA program for the specified duration
    simulation_config = SimulationConfig(duration=10_000)  # In clock cycles = 4ns
    job = qmm.simulate(config, cz, simulation_config)
    job.get_simulated_samples().con1.plot()
    plt.show()
else:
    # Open the quantum machine
    qm = qmm.open_qm(config)
    # Send the QUA program to the OPX, which compiles and executes it
    job = qm.execute(cz)
    # Prepare the figure for live plotting
    fig = plt.figure()
    interrupt_on_close(fig, job)
    # Tool to easily fetch results from the OPX (results_handle used in it)
    results = fetching_tool(job, ["n", "I1", "Q1", "I2", "Q2", "I3", "Q3", "I4", "Q4"], mode="live")
    # Live plotting
    while results.is_processing():
        # Fetch results
        n, I1, Q1, I2, Q2, I3, Q3, I4, Q4 = results.fetch_all()
        # Convert the results into Volts
        I1, Q1 = u.demod2volts(I1, readout_len), u.demod2volts(Q1, readout_len)
        I2, Q2 = u.demod2volts(I2, readout_len), u.demod2volts(Q2, readout_len)
        I3, Q3 = u.demod2volts(I3, readout_len), u.demod2volts(Q3, readout_len)
        I4, Q4 = u.demod2volts(I4, readout_len), u.demod2volts(Q4, readout_len)        
        # Progress bar
        progress_counter(n, n_avg, start_time=results.start_time)
        # Plot
        plt.suptitle(f"CZ chevron sweeping the flux on qubit {qubit_to_flux_tune}")
        # plt.subplot(221)
        # plt.cla()
        # plt.pcolor(amps * scale_reference + flux_offset, 4 * ts, I1)
        # plt.title("q1 - I [V]")
        # plt.ylabel("Interaction time (ns)")
        # plt.subplot(223)
        # plt.cla()
        # plt.pcolor(amps * scale_reference + flux_offset, 4 * ts, Q1)
        # plt.title("q1 - Q [V]")
        # plt.xlabel("Flux amplitude (V)")
        # plt.ylabel("Interaction time (ns)")
        # plt.subplot(222)
        # plt.cla()
        # plt.pcolor(amps * scale_reference + flux_offset, 4 * ts, I2)
        # plt.title("q2 - I [V]")
        # plt.subplot(224)
        # plt.cla()
        # plt.pcolor(amps * scale_reference + flux_offset, 4 * ts, Q2)
        # plt.title("q2 - Q [V]")
        # plt.xlabel("Flux amplitude (V)")

        plt.subplot(221)
        plt.cla()
        plt.pcolor(amps * scale_reference, 4 * ts, I3)
        plt.title("q3 - I [V]")
        plt.subplot(223)
        plt.cla()
        plt.pcolor(amps * scale_reference, 4 * ts, Q3)
        plt.title("q3 - Q [V]")
        plt.xlabel("Flux amplitude (V)")

        plt.subplot(222)
        plt.cla()
        plt.pcolor(amps * scale_reference, 4 * ts, I4)
        plt.title("q4 - I [V]")
        plt.subplot(224)
        plt.cla()
        plt.pcolor(amps * scale_reference, 4 * ts, Q4)
        plt.title("q4 - Q [V]")
        plt.xlabel("Flux amplitude (V)")        


        plt.tight_layout()
        plt.pause(0.1)
    # Close the quantum machines at the end in order to put all flux biases to 0 so that the fridge doesn't heat-up
    qm.close()
    plt.show()
    # np.savez(save_dir/'cz', I1=I1, Q1=Q1, I2=I2, Q2=Q2, ts=ts, amps=amps)