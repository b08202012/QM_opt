from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm import SimulationConfig
import matplotlib.pyplot as plt
from qualang_tools.loops import from_array
from qualang_tools.results import fetching_tool
from qualang_tools.plot import interrupt_on_close
from qualang_tools.results import progress_counter
from exp.RO_macros import multiRO_declare, multiRO_measurement, multiRO_pre_save
from qualang_tools.units import unit
import numpy as np
import time
import xarray as xr
import warnings
u = unit(coerce_to_integer=True)
warnings.filterwarnings("ignore")

def cz_gate(flux_Qi, flux_Ci,cz_amp,cz_time,c_amp):
    play("const" * amp(cz_amp),f"q{flux_Qi}_z", duration=cz_time/4)
    play("const" * amp(c_amp),f"q{flux_Ci}_z",duration=cz_time/4)

def cz_gate_compensate(flux_Qi, flux_Ci,cz_amp,cz_time,c_amp,control_Qi,target_Qi):
    play("const" * amp(cz_amp),f"q{flux_Qi}_z", duration=cz_time/4)
    play("const" * amp(c_amp),f"q{flux_Ci}_z",duration=cz_time/4)
    #frame_rotation_2pi(1.958/(2*np.pi), f"q{control_Qi}_xy")
    frame_rotation_2pi(-0.051673306854385946/(2*np.pi), f"q{target_Qi}_xy")

def CZ_phase_ramsey(cz_amps_range,cz_amps_resolution,cz_time,couplerz_amps_range,couplerz_amps_resolution,ro_element,flux_Qi,control_Qi,target_Qi,flux_Ci,preprocess,qmm,config,n_avg=1,initializer=None,simulate=True):
    cz_amps_array = np.arange(cz_amps_range[0],cz_amps_range[1],cz_amps_resolution)
    couplerz_amps_array = np.arange(couplerz_amps_range[0],couplerz_amps_range[1],couplerz_amps_resolution)
    phase_array = np.linspace(0, 1, 50)
    amps_len = len(cz_amps_array)
    camps_len = len(couplerz_amps_array)
    phase_array_len = len(phase_array)

    with program() as cz_phase:
        iqdata_stream = multiRO_declare( ro_element )
        n = declare(int)
        n_st = declare_stream()
        cz_amps = declare(fixed)
        couplerz_amps = declare(fixed)
        control = declare(int) # control qubit at 0 or 1
        phase = declare(fixed)  # second gate is x90 or y90 
        with for_(n, 0, n < n_avg, n + 1):
            with for_(*from_array(couplerz_amps, couplerz_amps_array)):
                with for_(*from_array(cz_amps, cz_amps_array)):
                    with for_each_( control, [0,1]):
                        with for_(*from_array(phase, phase_array)):
                            # initializaion
                            if initializer is None:
                                wait(1*u.us,ro_element)
                            else:
                                try:
                                    initializer[0](*initializer[1])
                                except:
                                    wait(1*u.us,ro_element)

                            # operation
                            with if_(control==1):
                                play("x180",f"q{control_Qi}_xy")
                            play("y90", f"q{target_Qi}_xy")
                            align(f"q{target_Qi}_xy",f"q{flux_Qi}_z",f'q{flux_Ci}_z')
                            wait(5)
                            cz_gate(flux_Qi,flux_Ci,cz_amps,cz_time,couplerz_amps)
                            align(f"q{flux_Qi}_z",f"q{target_Qi}_xy")
                            wait(5)
                            
                            frame_rotation_2pi(phase,f"q{target_Qi}_xy")
                            play("-y90", f"q{target_Qi}_xy")
                            align()
                            wait(5)
                                
                            # wait(cz_time*u.ns)
                            

                            # Readout
                            multiRO_measurement(iqdata_stream, ro_element, weights="rotated_")
            save(n, n_st)

        with stream_processing():
            n_st.save("iteration")
            match preprocess:
                case "shot":
                    multiRO_pre_save(iqdata_stream, ro_element, (n_avg,camps_len,amps_len,2,phase_array_len) ,stream_preprocess="shot")
                case _:
                    multiRO_pre_save(iqdata_stream, ro_element, (camps_len,amps_len,2,phase_array_len))
    if simulate:
        simulation_config = SimulationConfig(duration=20000)  # In clock cycles = 4ns
        job = qmm.simulate(config, cz_phase, simulation_config)
        job.get_simulated_samples().con1.plot()
        job.get_simulated_samples().con2.plot()
        plt.show()
    else:
        # Open the quantum machine
        qm = qmm.open_qm(config)
        # Send the QUA program to the OPX, which compiles and executes it
        job = qm.execute(cz_phase)
        # Get results from QUA program
        ro_ch_name = []
        for r_name in ro_element:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")
        data_list = ro_ch_name + ["iteration"]   
        results = fetching_tool(job, data_list=data_list, mode="live")

        # Live plotting
        while results.is_processing():
            # Fetch results
            fetch_data = results.fetch_all()
            # Progress bar
            iteration = fetch_data[-1]
            progress_counter(iteration, n_avg, start_time=results.start_time)
            # Plot
            plt.tight_layout()
            time.sleep(1)

        # Measurement finished 
        fetch_data = results.fetch_all()
        qm.close()
        output_data = {}

        match preprocess:
            case "shot":
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","shot","c_amp","cz_amp","control","phase"],
                                np.squeeze(np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) ))
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]), "shot":np.arange(n_avg),"c_amp":couplerz_amps_array, "cz_amp": cz_amps_array, "control":[0,1], "phase":phase_array }
                )
            case _:
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","c_amp","cz_amp","control","phase"],
                                np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) )
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]), "c_amp":couplerz_amps_array, "cz_amp": cz_amps_array, "control":[0,1], "phase":phase_array }
                )

        return dataset



def CZ_phase_diff(cz_amps_range,cz_amps_resolution,cz_time,couplerz_amps_range,couplerz_amps_resolution,ro_element,flux_Qi,control_Qi,target_Qi,flux_Ci,preprocess,qmm,config,n_avg=1,initializer=None,simulate=True):
    cz_amps_array = np.arange(cz_amps_range[0],cz_amps_range[1],cz_amps_resolution)
    couplerz_amps_array = np.arange(couplerz_amps_range[0],couplerz_amps_range[1],couplerz_amps_resolution)
    amps_len = len(cz_amps_array)
    camps_len = len(couplerz_amps_array)

    with program() as cz_phase:
        iqdata_stream = multiRO_declare( ro_element )
        n = declare(int)
        n_st = declare_stream()
        cz_amps = declare(fixed)
        couplerz_amps = declare(fixed)
        control = declare(int) # control qubit at 0 or 1
        rotate = declare(int)  # second gate is x90 or y90 
        with for_(n, 0, n < n_avg, n + 1):
            with for_(*from_array(couplerz_amps, couplerz_amps_array)):
                with for_(*from_array(cz_amps, cz_amps_array)):
                    with for_each_( control, [0,1]):

                        with for_each_(rotate, [0,1, 2]):

                            # initializaion
                            if initializer is None:
                                wait(1*u.us,ro_element)
                            else:
                                try:
                                    initializer[0](*initializer[1])
                                except:
                                    wait(1*u.us,ro_element)

                            # operation
                            with if_(control==1):
                                play("x180",f"q{control_Qi}_xy")

                            play("y90", f"q{target_Qi}_xy")

                            align(f"q{target_Qi}_xy",f"q{flux_Qi}_z",f'q{flux_Ci}_z')
                            wait(20*u.ns)
                            cz_gate(flux_Qi,flux_Ci,cz_amps,cz_time,couplerz_amps)
                            align(f"q{flux_Qi}_z",f"q{target_Qi}_xy")
                            wait(20*u.ns)
                            with if_(rotate==1):
                                play("x90", f"q{target_Qi}_xy")
                                align()
                                wait(20*u.ns)
                            with elif_(rotate == 0):
                                play("-y90", f"q{target_Qi}_xy")
                                align()
                                wait(20*u.ns)
                                
                            # wait(cz_time*u.ns)
                            


                            # Readout
                            multiRO_measurement(iqdata_stream, ro_element, weights="rotated_")
            save(n, n_st)

        with stream_processing():
            n_st.save("iteration")
            match preprocess:
                case "shot":

                    multiRO_pre_save(iqdata_stream, ro_element, (n_avg,camps_len,amps_len,2,3) ,stream_preprocess="shot")
                case _:
                    multiRO_pre_save(iqdata_stream, ro_element, (camps_len,amps_len,2,3))

    if simulate:
        simulation_config = SimulationConfig(duration=20000)  # In clock cycles = 4ns
        job = qmm.simulate(config, cz_phase, simulation_config)
        job.get_simulated_samples().con1.plot()
        job.get_simulated_samples().con2.plot()
        plt.show()
    else:
        # Open the quantum machine
        qm = qmm.open_qm(config)
        # Send the QUA program to the OPX, which compiles and executes it
        job = qm.execute(cz_phase)
        # Get results from QUA program
        ro_ch_name = []
        for r_name in ro_element:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")
        data_list = ro_ch_name + ["iteration"]   
        results = fetching_tool(job, data_list=data_list, mode="live")

        # Live plotting
        while results.is_processing():
            # Fetch results
            fetch_data = results.fetch_all()
            # Progress bar
            iteration = fetch_data[-1]
            progress_counter(iteration, n_avg, start_time=results.start_time)
            # Plot
            plt.tight_layout()
            time.sleep(1)

        # Measurement finished 
        fetch_data = results.fetch_all()
        qm.close()
        output_data = {}

        match preprocess:
            case "shot":
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","shot","c_amp","cz_amp","control","rotate"],
                                np.squeeze(np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) ))
                dataset = xr.Dataset(
                    output_data,

                    coords={ "mixer":np.array(["I","Q"]), "shot":np.arange(n_avg),"c_amp":couplerz_amps_array, "cz_amp": cz_amps_array, "control":[0,1], "rotate":[0,1, 2] }

                )
            case _:
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","c_amp","cz_amp","control","rotate"],
                                np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) )
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]), "c_amp":couplerz_amps_array, "cz_amp": cz_amps_array, "control":[0,1], "rotate":[0,1, 2] }

                )

        return dataset
    
def CZ_phase_diff_1D(cz_amps,cz_time,ro_element,flux_Qi,control_Qi,target_Qi,flux_Ci,preprocess,qmm,config,n_avg=1,initializer=None,simulate=True):

    with program() as cz_phase:
        iqdata_stream = multiRO_declare( ro_element )
        n = declare(int)
        n_st = declare_stream()
        control = declare(int) # control qubit at 0 or 1
        rotate = declare(int)  # second gate is x90 or y90 
        with for_(n, 0, n < n_avg, n + 1):
            with for_each_( control, [0,1]):
                with for_each_(rotate, [0,1]):
                    # initializaion
                    if initializer is None:
                        wait(1*u.us,ro_element)
                    else:
                        try:
                            initializer[0](*initializer[1])
                        except:
                            wait(1*u.us,ro_element)

                    # operation
                    with if_(control==1):
                        play("x180",f"q{control_Qi}_xy")
                    play("x90", f"q{target_Qi}_xy")
                    align(f"q{target_Qi}_xy",f"q{flux_Qi}_z",f'q{flux_Ci}_z')
                    wait(20*u.ns)
                    cz_gate(flux_Qi,flux_Ci,cz_amps,cz_time,0.0)
                    align(f"q{flux_Qi}_z",f"q{target_Qi}_xy")
                    wait(20*u.ns)
                    with if_(rotate==1):
                        play("-y90", f"q{target_Qi}_xy")
                    with else_():
                        play("x90", f"q{target_Qi}_xy")
                    wait(cz_time*u.ns)
                    wait(100*u.ns)

                    # Readout
                    multiRO_measurement(iqdata_stream, ro_element, weights="rotated_")
            save(n, n_st)

        with stream_processing():
            n_st.save("iteration")
            match preprocess:
                case "shot":
                    multiRO_pre_save(iqdata_stream, ro_element, (n_avg,2,2) ,stream_preprocess="shot")
                case _:
                    multiRO_pre_save(iqdata_stream, ro_element, (2,2))
    if simulate:
        simulation_config = SimulationConfig(duration=20000)  # In clock cycles = 4ns
        job = qmm.simulate(config, cz_phase, simulation_config)
        job.get_simulated_samples().con1.plot()
        job.get_simulated_samples().con2.plot()
        plt.show()
    else:
        # Open the quantum machine
        qm = qmm.open_qm(config)
        # Send the QUA program to the OPX, which compiles and executes it
        job = qm.execute(cz_phase)
        # Get results from QUA program
        ro_ch_name = []
        for r_name in ro_element:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")
        data_list = ro_ch_name + ["iteration"]   
        results = fetching_tool(job, data_list=data_list, mode="live")

        # Live plotting
        while results.is_processing():
            # Fetch results
            fetch_data = results.fetch_all()
            # Progress bar
            iteration = fetch_data[-1]
            progress_counter(iteration, n_avg, start_time=results.start_time)
            # Plot
            plt.tight_layout()
            time.sleep(1)

        # Measurement finished 
        fetch_data = results.fetch_all()
        qm.close()
        output_data = {}

        match preprocess:
            case "shot":
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","shot","control","rotate"],
                                np.squeeze(np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) ))
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]), "shot":np.arange(n_avg), "control":[0,1], "rotate":[0,1] }
                )
            case _:
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","control","rotate"],
                                np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) )
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]), "control":[0,1], "rotate":[0,1] }
                )
        return dataset

def CZ_phase_diff_time(cz_amps_range,cz_amps_resolution,cz_time_max,cz_time_resolution,c_amps,ro_element,flux_Qi,control_Qi,target_Qi,flux_Ci,qmm,config,n_avg=1,initializer=None,simulate=True):
    cz_amps_array = np.arange(cz_amps_range[0],cz_amps_range[1],cz_amps_resolution)
    cc_resolution = (cz_time_resolution/4.)*u.ns
    cc_max_qua = (cz_time_max/4.)*u.ns
    cc_qua = np.arange( 4, cc_max_qua, cc_resolution)
    print(cc_qua)
    cz_time_array = cc_qua*4

    amps_len = len(cz_amps_array)
    time_len = len(cz_time_array)

    with program() as cz_phase:
        iqdata_stream = multiRO_declare( ro_element )
        n = declare(int)
        n_st = declare_stream()
        cz_amps = declare(fixed)
        cc = declare(int)
        control = declare(int) # control qubit at 0 or 1
        rotate = declare(int)  # second gate is x90 or y90 
        with for_(n, 0, n < n_avg, n + 1):
            with for_(*from_array(cz_amps, cz_amps_array)):
                with for_(*from_array(cc, cc_qua)):
                    with for_each_( control, [0,1]):
                        with for_each_(rotate, [0,1]):
                            # initializaion
                            if initializer is None:
                                wait(1*u.us,ro_element)
                            else:
                                try:
                                    initializer[0](*initializer[1])
                                except:
                                    wait(1*u.us,ro_element)

                            # operation
                            with if_(control==1):
                                play("x180",f"q{control_Qi}_xy")
                            play("x90", f"q{target_Qi}_xy")
                            align(f"q{target_Qi}_xy",f"q{flux_Qi}_z",f'q{flux_Ci}_z')
                            wait(40*u.ns)
                            cz_gate(flux_Qi,flux_Ci,cz_amps,cc*4,c_amps)
                            #cz_gate_compensate(flux_Qi, flux_Ci,cz_amps,cc*4,c_amps,control_Qi,target_Qi)
                            with if_(rotate==1):
                                frame_rotation_2pi(0.25, f"q{target_Qi}_xy")
                            align(f"q{flux_Qi}_z",f"q{target_Qi}_xy")
                            wait(20 *u.ns)
                            play("x90", f"q{target_Qi}_xy")
                            wait(cc)
                            wait(140*u.ns)

                            # Readout
                            multiRO_measurement(iqdata_stream, ro_element, weights="rotated_")
            save(n, n_st)

        with stream_processing():
            n_st.save("iteration")
            multiRO_pre_save(iqdata_stream, ro_element, (amps_len,time_len,2,2) )
    if simulate:
        simulation_config = SimulationConfig(duration=20000)  # In clock cycles = 4ns
        job = qmm.simulate(config, cz_phase, simulation_config)
        job.get_simulated_samples().con1.plot()
        job.get_simulated_samples().con2.plot()
        plt.show()
    else:
        # Open the quantum machine
        qm = qmm.open_qm(config)
        # Send the QUA program to the OPX, which compiles and executes it
        job = qm.execute(cz_phase)
        # Get results from QUA program
        ro_ch_name = []
        for r_name in ro_element:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")
        data_list = ro_ch_name + ["iteration"]   
        results = fetching_tool(job, data_list=data_list, mode="live")

        # Live plotting
        while results.is_processing():
            # Fetch results
            fetch_data = results.fetch_all()
            # Progress bar
            iteration = fetch_data[-1]
            progress_counter(iteration, n_avg, start_time=results.start_time)
            # Plot
            plt.tight_layout()
            time.sleep(1)

        # Measurement finished 
        fetch_data = results.fetch_all()
        qm.close()
        output_data = {}

        for r_idx, r_name in enumerate(ro_element):
            output_data[r_name] = ( ["mixer","cz_amp","cz_time","control","rotate"],
                                np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) )
        dataset = xr.Dataset(
            output_data,
            coords={ "mixer":np.array(["I","Q"]),"cz_amp":cz_amps_array,"cz_time":cz_time_array,"control":np.array([0,1]), "rotate": np.array([0,1])}
        )

        return dataset

def CZ_phase_compensate(c_amp,cz_amp,cz_time,ro_element,flux_Qi,target_Qi,flux_Ci,preprocess,qmm,config,n_avg=1,initializer=None,simulate=True):
    with program() as cz_phase:
        iqdata_stream = multiRO_declare( ro_element )
        n = declare(int)
        n_st = declare_stream()
        rotate = declare(int)  # second gate is x90 or y90 
        with for_(n, 0, n < n_avg, n + 1):
            with for_each_(rotate, [0,1]):
                # initializaion
                if initializer is None:
                    wait(1*u.us,ro_element)
                else:
                    try:
                        initializer[0](*initializer[1])
                    except:
                        wait(1*u.us,ro_element)

                # operation
                play("x90", f"q{target_Qi}_xy")
                align(f"q{target_Qi}_xy",f"q{flux_Qi}_z",f"q{flux_Ci}_z")
                wait(20*u.ns)
                cz_gate(flux_Qi,flux_Ci,cz_amp,cz_time,c_amp)

                align(f"q{flux_Qi}_z",f"q{target_Qi}_xy")
                wait(20*u.ns)
                with if_(rotate==1):
                    play("-y90", f"q{target_Qi}_xy")
                with else_():
                    play("x90", f"q{target_Qi}_xy")
                wait(cz_time*u.ns)
                wait(100*u.ns)

                # Readout
                multiRO_measurement(iqdata_stream, ro_element, weights="rotated_")
            save(n, n_st)

        with stream_processing():
            n_st.save("iteration")
            match preprocess:
                case "shot":
                    multiRO_pre_save(iqdata_stream, ro_element, (n_avg,2) ,stream_preprocess="shot")
                case _:
                    multiRO_pre_save(iqdata_stream, ro_element, (2,))
    if simulate:
        simulation_config = SimulationConfig(duration=20_000)  # In clock cycles = 4ns
        job = qmm.simulate(config, cz_phase, simulation_config)
        job.get_simulated_samples().con1.plot()
        # job.get_simulated_samples().con2.plot()
        plt.show()
    else:
        # Open the quantum machine
        qm = qmm.open_qm(config)
        # Send the QUA program to the OPX, which compiles and executes it
        job = qm.execute(cz_phase)
        # Get results from QUA program
        ro_ch_name = []
        for r_name in ro_element:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")
        data_list = ro_ch_name + ["iteration"]   
        results = fetching_tool(job, data_list=data_list, mode="live")

        # Live plotting
        while results.is_processing():
            # Fetch results
            fetch_data = results.fetch_all()
            # Progress bar
            iteration = fetch_data[-1]
            progress_counter(iteration, n_avg, start_time=results.start_time)
            # Plot
            plt.tight_layout()
            time.sleep(1)

        # Measurement finished 
        fetch_data = results.fetch_all()
        qm.close()
        output_data = {}

        match preprocess:
            case "shot":
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","shot","rotate"],
                                np.squeeze(np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) ))
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]), "shot":np.arange(n_avg), "rotate":[0,1] }
                )
            case _:
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","rotate"],
                                np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) )
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]),"rotate":[0,1] }
                )

        return dataset

def CZ_fidelity(control_q_phase, target_q_phase, c_amp,cz_amp,cz_time,ro_element,flux_Qi, control_Qi,target_Qi,flux_Ci,preprocess,qmm,config,n_avg=1,initializer=None,simulate=True):
    with program() as cz_phase:
        iqdata_stream = multiRO_declare( ro_element )
        n = declare(int)
        n_st = declare_stream()
        rotate = declare(int)  # second gate is x90 or y90 
        with for_(n, 0, n < n_avg, n + 1):
            # initializaion
            if initializer is None:
                wait(1*u.us,ro_element)
            else:
                try:
                    initializer[0](*initializer[1])
                except:
                    wait(1*u.us,ro_element)

            # operation
            play("x180", f"q{target_Qi}_xy")
            play("x180", f"q{control_Qi}_xy")

            play("y90", f"q{target_Qi}_xy")
            play("y90", f"q{control_Qi}_xy")
            align()
            wait(5)     

            # wait(25)
            play("const" * amp(cz_amp),f"q{flux_Qi}_z")
            play("const" * amp(c_amp),f"q{flux_Ci}_z")
            align()
            frame_rotation_2pi(control_q_phase/(2*np.pi), f"q{control_Qi}_xy")
            frame_rotation_2pi(target_q_phase/(2*np.pi), f"q{target_Qi}_xy")

            align()
            wait(5)
            play("x180", f"q{target_Qi}_xy")
            play("y90", f"q{target_Qi}_xy")
            # play("x90", f"q{control_Qi}_xy")
            align()
            # Readout
            multiRO_measurement(iqdata_stream, ro_element, weights="rotated_")
            save(n, n_st)

        with stream_processing():
            n_st.save("iteration")
            match preprocess:
                case "shot":
                    multiRO_pre_save(iqdata_stream, ro_element, (n_avg,) ,stream_preprocess="shot")
                case _:
                    multiRO_pre_save(iqdata_stream, ro_element, ())
    if simulate:
        simulation_config = SimulationConfig(duration=20_000)  # In clock cycles = 4ns
        job = qmm.simulate(config, cz_phase, simulation_config)
        job.get_simulated_samples().con1.plot()
        # job.get_simulated_samples().con2.plot()

        plt.show()
    else:
        # Open the quantum machine
        qm = qmm.open_qm(config)
        # Send the QUA program to the OPX, which compiles and executes it
        job = qm.execute(cz_phase)
        # Get results from QUA program
        ro_ch_name = []
        for r_name in ro_element:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")
        data_list = ro_ch_name + ["iteration"]   
        results = fetching_tool(job, data_list=data_list, mode="live")

        # Live plotting
        while results.is_processing():
            # Fetch results
            fetch_data = results.fetch_all()
            # Progress bar
            iteration = fetch_data[-1]
            progress_counter(iteration, n_avg, start_time=results.start_time)
            # Plot
            plt.tight_layout()
            time.sleep(1)

        # Measurement finished 
        fetch_data = results.fetch_all()
        qm.close()
        output_data = {}

        match preprocess:
            case "shot":
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer","shot"],
                                np.squeeze(np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) ))
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]), "shot":np.arange(n_avg) }
                )
            case _:
                for r_idx, r_name in enumerate(ro_element):
                    output_data[r_name] = ( ["mixer"],
                                np.array([fetch_data[r_idx*2], fetch_data[r_idx*2+1]]) )
                dataset = xr.Dataset(
                    output_data,
                    coords={ "mixer":np.array(["I","Q"]) }
                )

        return dataset


def plot_cz_phase_diff(x,y,ax=None):
    """
    x in (N,)
    y in (2,N)
    2 is control qubit at 0 or 1
    N is second x90 gate rotated angle
    """
    if ax == None:
        fig, ax = plt.subplots()
    ax.set_title("CZ phase difference")
    ax.set_xlabel("angle (rad)")
    ax.plot(x,y[0],"-",color='b',label="control at 0")
    ax.plot(x,y[1],"--",color="r",label="control at 1")
