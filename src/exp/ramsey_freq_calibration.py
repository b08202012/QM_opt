

from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm import SimulationConfig
# from configuration import *
import matplotlib.pyplot as plt
from qualang_tools.loops import from_array
from qualang_tools.results import fetching_tool, progress_counter
from qualang_tools.plot import interrupt_on_close
from exp.RO_macros import multiRO_declare, multiRO_measurement, multiRO_pre_save
from qualang_tools.plot.fitting import Fit
import xarray as xr
import warnings
warnings.filterwarnings("ignore")
from qualang_tools.units import unit


#######################
# AUXILIARY FUNCTIONS #
#######################
u = unit(coerce_to_integer=True)

###################
# The QUA program #
###################

def ramsey_freq_calibration( virtial_detune_freq, q_name:list, ro_element:list, config, qmm:QuantumMachinesManager, n_avg:int=100, simulate = False, initializer:tuple=None ):
    """
    Use positive and nagative detuning refence to freq in config to get measured ramsey oscillation frequency.
    evo_time unit is tick (4ns)
    virtial_detune_freq unit in MHz can't larger than 2
    """
    point_per_period = 20
    Ramsey_period = (1e3/virtial_detune_freq)* u.ns
    tick_resolution = (Ramsey_period//(4*point_per_period))
    evo_time_tick_max = tick_resolution *point_per_period*6
    print(f"time resolution {tick_resolution*4} ,max time {evo_time_tick_max*4}")
    evo_time_tick = np.arange( 4, evo_time_tick_max, tick_resolution)
    evo_time = evo_time_tick*4
    time_len = len(evo_time)
    with program() as ramsey:
        iqdata_stream = multiRO_declare( ro_element )
        n = declare(int)
        n_st = declare_stream()
        t = declare(int)  # QUA variable for the idle time
        phi = declare(fixed)  # Phase to apply the virtual Z-rotation
        phi_idx = declare(bool,)
        with for_(n, 0, n < n_avg, n + 1):
            with for_each_( phi_idx, [True, False]):
                with for_(*from_array(t, evo_time_tick)):

                    # Rotate the frame of the second x90 gate to implement a virtual Z-rotation
                    # 4*tau because tau was in clock cycles and 1e-9 because tau is ns
                    
                    # Init
                    if initializer is None:
                        wait(100*u.us)
                        #wait(thermalization_time * u.ns)
                    else:
                        try:
                            initializer[0](*initializer[1])
                        except:
                            print("Initializer didn't work!")
                            wait(100*u.us)

                    # Operation
                    True_value = Cast.mul_fixed_by_int(virtial_detune_freq * 1e-3, 4 * t)
                    False_value = Cast.mul_fixed_by_int(-virtial_detune_freq * 1e-3, 4 * t)
                    assign(phi, Util.cond(phi_idx, True_value, False_value))

                    for q in q_name:
                        play("x90", q)  # 1st x90 gate

                    for q in q_name:
                        wait(t, q)

                    for q in q_name:
                        frame_rotation_2pi(phi, q)  # Virtual Z-rotation
                        play("x90", q)  # 2st x90 gate

                    # Align after playing the qubit pulses.
                    align()
                    # Readout
                    multiRO_measurement(iqdata_stream, ro_element, weights="rotated_")         
                

            # Save the averaging iteration to get the progress bar
            save(n, n_st)

        with stream_processing():
            n_st.save("iteration")
            multiRO_pre_save(iqdata_stream, ro_element, (2,time_len) )

    ###########################
    # Run or Simulate Program #
    ###########################


    if simulate:
        # Simulates the QUA program for the specified duration
        simulation_config = SimulationConfig(duration=20_000)  # In clock cycles = 4ns
        job = qmm.simulate(config, ramsey, simulation_config)
        job.get_simulated_samples().con1.plot()
        job.get_simulated_samples().con2.plot()
        plt.show()

    else:
        # Open the quantum machine
        qm = qmm.open_qm(config)
        # Send the QUA program to the OPX, which compiles and executes it
        job = qm.execute(ramsey)
        # Get results from QUA program
        ro_ch_name = []
        for r_name in ro_element:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")
        data_list = ro_ch_name + ["iteration"]   
        results = fetching_tool(job, data_list=data_list, mode="live")
        # Live plotting

        fig, ax = plt.subplots(2, len(ro_element))
        interrupt_on_close(fig, job)  # Interrupts the job when closing the figure
        fig.suptitle("Frequency calibration")

        # Live plotting
        while results.is_processing():
            # Fetch results
            fetch_data = results.fetch_all()

            # Progress bar
            iteration = fetch_data[-1]
            progress_counter(iteration, n_avg, start_time=results.start_time)

        # Close the quantum machines at the end in order to put all flux biases to 0 so that the fridge doesn't heat-up
        qm.close()
        # Creating an xarray dataset
        output_data = {}
        for r_idx, r_name in enumerate(ro_element):
            output_data[r_name] = ( ["mixer","frequency","time"],
                                    np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) )
        dataset = xr.Dataset(
            output_data,
            coords={ "mixer":np.array(["I","Q"]), "frequency": np.array([virtial_detune_freq,-virtial_detune_freq]), "time":evo_time }
        )
        # dataset.attrs["ref_xy_IF"] = ref_xy_IF
        # dataset.attrs["ref_xy_LO"] = ref_xy_LO

        return dataset
        
def plot_dual_Ramsey_oscillation( x, y, ax=None ):
    """
    y in shape (2,N)
    2 is postive and negative
    N is evo_time_point
    """
    if ax == None:
        fig, ax = plt.subplots()
    ax.plot(x, y[0], "o",label="positive")
    ax.plot(x, y[1], "o",label="negative")
    ax.set_xlabel("Free Evolution Times [ns]")
    ax.legend()

    if ax == None:
        return fig
    
def plot_ana_result( evo_time, data, detuning, ax=None ):
    """
    data in shape (2,N)
    2 is postive and negative
    N is evo_time_point
    """
    if ax == None:
        fig, ax = plt.subplots()
    fit = Fit()
    plot_dual_Ramsey_oscillation(evo_time, data, ax)
    ax.set_title(f"ZZ-Ramsey measurement with virtual detuning {detuning} MHz")

    ana_dict_pos = fit.ramsey(evo_time, data[0], plot=False)
    ana_dict_neg = fit.ramsey(evo_time, data[1], plot=False)

    ax.set_xlabel("Idle times [ns]")

    freq_pos = ana_dict_pos['f'][0]*1e3
    freq_neg = ana_dict_neg['f'][0]*1e3
    ax.plot(evo_time, ana_dict_pos["fit_func"](evo_time), label=f"Positive freq: {freq_pos:.3f} MHz")
    ax.plot(evo_time, ana_dict_neg["fit_func"](evo_time), label=f"Negative freq: {freq_neg:.3f} MHz")
    ax.text(0.07, 0.9, f"Real Detuning freq : {(freq_pos-freq_neg)/2:.3f}", fontsize=10, transform=ax.transAxes)

    ax.legend()
    plt.tight_layout()
    return (freq_pos-freq_neg)/2

if __name__ == '__main__':
    qmm = QuantumMachinesManager(host=qop_ip, port=qop_port, cluster_name=cluster_name, octave=octave_config)
    n_avg = 1000  # Number of averages


    ro_element = ["rr1"]
    q_name =  ["q1_xy"]
    virtual_detune = 1 # Unit in MHz
    output_data, evo_time = Ramsey_freq_calibration( virtual_detune, q_name, ro_element, config, qmm, n_avg=n_avg, simulate=False)
    #   Data Saving   # 
    save_data = False
    if save_data:
        from save_data import save_npz
        import sys
        save_progam_name = sys.argv[0].split('\\')[-1].split('.')[0]  # get the name of current running .py program
        save_npz(save_dir, save_progam_name, output_data)

    plot_ana_result(evo_time,output_data[ro_element[0]][0],virtual_detune)
    # # Plot
    plt.show()

