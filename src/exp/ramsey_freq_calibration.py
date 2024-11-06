

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
u = unit(coerce_to_integer=True)

import exp.config_par as gc
import xarray as xr
from exp.QMMeasurement import QMMeasurement

class RamseyFreqCalibration( QMMeasurement ):
    """
    Use positive and nagative detuning refence to freq in config to get measured ramsey oscillation frequency.
    evo_time unit is tick (4ns)
    virtial_detune_freq unit in MHz can't larger than 2

    return: \n
    dataset \n
    coors: ["mixer","frequency","time"]\n
    attrs: ref_xy_IF, ref_xy_LO, z_offset\n
    """

    def __init__( self, config, qmm: QuantumMachinesManager ):
        super().__init__( config, qmm )

        self.ro_elements = ["q4_ro"]
        self.z_elements = ["q1_z"]
        self.xy_elements = ["q4_xy"]
        
        self.preprocess = "ave"
        self.initializer = None
        
        # self.sweep_type = "z_pulse"


        self.virtial_detune_freq = 5

        self.point_per_period = 20
        self.max_period = 6
        

    def _get_qua_program( self ):
        

        self.qua_cc_evo = self._lin_time_array()
        print(self.qua_cc_evo)
        self._attribute_config()

        with program() as ramsey:
            iqdata_stream = multiRO_declare( self.ro_elements )
            n = declare(int)
            n_st = declare_stream()
            t = declare(int)  # QUA variable for the idle time
            phi = declare(fixed)  # Phase to apply the virtual Z-rotation
            phi_idx = declare(bool,)
            with for_(n, 0, n < self.shot_num, n + 1):
                with for_each_( phi_idx, [True, False]):
                    with for_(*from_array(t, self.qua_cc_evo)):

                        # Rotate the frame of the second x90 gate to implement a virtual Z-rotation
                        # 4*tau because tau was in clock cycles and 1e-9 because tau is ns
                        
                        # Init
                        if self.initializer is None:
                            wait(100*u.us)
                            #wait(thermalization_time * u.ns)
                        else:
                            try:
                                self.initializer[0](*self.initializer[1])
                            except:
                                print("Initializer didn't work!")
                                wait(100*u.us)

                        # Operation
                        True_value = Cast.mul_fixed_by_int(self.virtial_detune_freq * 1e-3, 4 * t)
                        False_value = Cast.mul_fixed_by_int(-self.virtial_detune_freq * 1e-3, 4 * t)
                        assign(phi, Util.cond(phi_idx, True_value, False_value))

                        for q in self.xy_elements:
                            play("x90", q)  # 1st x90 gate

                        for q in self.xy_elements:
                            wait(t, q)
                            #play("x180",q)
                            #wait(t,q)

                        for q in self.xy_elements:
                            frame_rotation_2pi(phi, q)  # Virtual Z-rotation
                            play("x90", q)  # 2st x90 gate

                        # Align after playing the qubit pulses.
                        wait(t)
                        wait(120*u.ns)
                        #align()
                        # Readout
                        multiRO_measurement(iqdata_stream, self.ro_elements, weights="rotated_")         
                    

                # Save the averaging iteration to get the progress bar
                save(n, n_st)

            with stream_processing():
                n_st.save("iteration")
                multiRO_pre_save(iqdata_stream, self.ro_elements, (2,len(self.qua_cc_evo)) )
        
        return ramsey
        
    
    def _get_fetch_data_list( self ):
        ro_ch_name = []
        for r_name in self.ro_elements:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")

        data_list = ro_ch_name + ["iteration"]   
        return data_list
    
    def _data_formation( self ):
        time_ns = self.qua_cc_evo*4

        coords = { 
            "mixer":np.array(["I","Q"]), 
            "frequency": np.array([self.virtial_detune_freq,-self.virtial_detune_freq]), 
            "time": time_ns
            }
        match self.preprocess:
            case "shot":
                dims_order = ["mixer","frequency","time"]
                coords["shot"] = np.arange(self.shot_num)
            case _:
                dims_order = ["mixer","frequency","time"]

        output_data = {}
        for r_idx, r_name in enumerate(self.ro_elements):
            data_array = np.array([ self.fetch_data[r_idx*2], self.fetch_data[r_idx*2+1]])
            output_data[r_name] = ( dims_order, np.squeeze(data_array))

        dataset = xr.Dataset( output_data, coords=coords )

        # dataset = dataset.transpose("mixer", "prepare_state", "frequency", "amp_ratio")

        self._attribute_config()
        dataset.attrs["ro_LO"] = self.ref_ro_LO
        dataset.attrs["ro_IF"] = self.ref_ro_IF
        dataset.attrs["xy_LO"] = self.ref_xy_LO
        dataset.attrs["xy_IF"] = self.ref_xy_IF
        dataset.attrs["z_offset"] = self.z_offset

        return dataset

    def _attribute_config( self ):
        self.ref_ro_IF = []
        self.ref_ro_LO = []
        for r in self.ro_elements:
            self.ref_ro_IF.append(gc.get_IF(r, self.config))
            self.ref_ro_LO.append(gc.get_LO(r, self.config))

        self.ref_xy_IF = []
        self.ref_xy_LO = []
        for xy in self.xy_elements:
            self.ref_xy_IF.append(gc.get_IF(xy, self.config))
            self.ref_xy_LO.append(gc.get_LO(xy, self.config))

        self.z_offset = []
        for z in self.z_elements:
            self.z_offset.append( gc.get_offset(z, self.config ))

    def _lin_time_array( self ):

        ramsey_period = abs(1e3/ self.virtial_detune_freq )* u.ns
        qua_cc_resolution = (ramsey_period//(4*self.point_per_period))
        if qua_cc_resolution < 1:
            print("Warning qua_cc_resolution <1 force to 1, virtial_detune_freq is to large or point_per_period is too large.")
            qua_cc_resolution = 1

        qua_cc_max_evo = qua_cc_resolution *self.point_per_period* self.max_period
        # print(f"time resolution {qua_cc_resolution*4} ,max time {evo_time_tick_max*4}")
        qua_cc_evo = np.arange( 4, qua_cc_max_evo, qua_cc_resolution)
        # evo_time = evo_time_tick*4
        
        return qua_cc_evo



        
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

