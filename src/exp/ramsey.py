from qm.QuantumMachinesManager import QuantumMachinesManager
from qm.qua import *
from qm import SimulationConfig
import matplotlib.pyplot as plt
from qualang_tools.loops import from_array
from qualang_tools.results import fetching_tool, progress_counter
from qualang_tools.plot import interrupt_on_close
from exp.RO_macros import multiRO_declare, multiRO_measurement, multiRO_pre_save
from qualang_tools.plot.fitting import Fit
# from common_fitting_func import *
import warnings
warnings.filterwarnings("ignore")
from qualang_tools.units import unit
u = unit(coerce_to_integer=True)
import xarray as xr
import time
from exp.QMMeasurement import QMMeasurement

class Ramsey( QMMeasurement ):
    def __init__( self, config, qmm: QuantumMachinesManager):
        super().__init__( config, qmm )
        """
        virtual_detune unit in MHz.\n
        time_max unit in us.\n
        time_resolution unit in us.\n
        """

        self.ro_elements = ["q0_ro"]
        self.xy_elements = ["q0_xy"]
        self.virtual_detune = 0
        self.max_time = 5
        self.time_resolution = 10
        self.n_avg = 100
        self.initializer = None
        self.simulate = False

    def _get_qua_program( self ):
        v_detune_qua = self.virtual_detune *u.MHz
        print(self.virtual_detune)
        cc_resolution = (self.time_resolution/4.) *u.us
        cc_max_qua = (self.max_time/4.) *u.us
        cc_qua = np.arange( 4, cc_max_qua, cc_resolution)
        print(cc_qua)

        self.evo_time = cc_qua*4
        time_len = len(cc_qua)
        with program() as ramsey:
            iqdata_stream = multiRO_declare( self.ro_elements )
            n = declare(int)
            n_st = declare_stream()
            cc = declare(int)  # QUA variable for the idle time, unit in clock cycle
            phi = declare(fixed)  # Phase to apply the virtual Z-rotation
            with for_(n, 0, n < self.n_avg, n + 1):
                with for_( *from_array(cc, cc_qua) ):
                    
                        # Init
                        if self.initializer is None:
                            wait(100*u.us)
                        else:
                            try:
                                self.initializer[0](*self.initializer[1])
                            except:
                                print("Initializer didn't work!")
                                wait(100*u.us)

                        # Operation
                        phi = Cast.mul_fixed_by_int( self.virtual_detune/1e3, 4 *cc)
                        # True_value =  v_detune_qua*4*cc
                        # False_value = v_detune_qua*4*cc

                        for xy in self.xy_elements:
                            play("x90", xy)  # 1st x90 gate
                            wait(cc, xy)
                            frame_rotation_2pi(phi, xy)  # Virtual Z-rotation
                            play("x90", xy)  # 2st x90 gate

                        # Align after playing the qubit pulses.
                        align()
                        # Readout
                        multiRO_measurement(iqdata_stream, self.ro_elements, weights="rotated_")         
                    

                # Save the averaging iteration to get the progress bar
                save(n, n_st)

            with stream_processing():
                n_st.save("iteration")
                multiRO_pre_save(iqdata_stream, self.ro_elements, (time_len,) )
        return ramsey
    
    def _get_fetch_data_list( self ):
        ro_ch_name = []
        for r_name in self.ro_elements:
            ro_ch_name.append(f"{r_name}_I")
            ro_ch_name.append(f"{r_name}_Q")

        data_list = ro_ch_name + ["iteration"]   
        return data_list
    
    def _data_formation( self ):
        output_data = {}

        for r_idx, r_name in enumerate(self.ro_elements):
            output_data[r_name] = ( ["mixer","time"],
                                np.array([self.fetch_data[r_idx*2], self.fetch_data[r_idx*2+1]]) )
        dataset = xr.Dataset(
            output_data,
            coords={ "mixer":np.array(["I","Q"]), "time": self.evo_time }
        )

        return dataset
       
def T2_fitting(time,signal):
    try:
        fit = Fit()
        decay_fit = fit.ramsey(time, signal, plot=False)
        qubit_T2 = np.round(np.abs(decay_fit["T2"][0]) / 4) * 4
    except Exception as e:     
        print(f"An error occurred: {e}")  
        qubit_T2 = 0
    return qubit_T2

def multi_T2_exp( repeat, time_max,time_resolution,ro_element,xy_elements,n_avg,config,qmm,virtual_detune=0,initializer=None ):

    raw_data = {}
    repetition = np.arange(repeat)
    for r in ro_element:
        raw_data[r] = []

    for i in range(repeat):
        print(f"{i}th T2")
        dataset = exp_ramsey( time_max,time_resolution,ro_element,xy_elements,n_avg,config,qmm,virtual_detune=virtual_detune,initializer=initializer )
        time = dataset.coords["time"].values
        for ro_name, data in dataset.data_vars.items():
            raw_data[ro_name].append(data)

    output_data = {}
    for r in ro_element:
        output_data[r] = (["repetition","mixer","time"], np.array(raw_data[r]))

    dataset = xr.Dataset(
        output_data,
        coords={ "mixer":np.array(["I","Q"]), "time": time, "repetition": repetition }
    )    
    dataset = dataset.transpose("mixer","repetition","time")

    return dataset

def plot_ramsey_oscillation( x, y, ax=None ):
    """
    y in shape (2,N)
    2 is postive and negative
    N is evo_time_point
    """
    if ax == None:
        fig, ax = plt.subplots()
    ax.plot(x, y, "-",label="T2")
    T2 = T2_fitting(x,y)
    print(T2)
    ax.set_xlabel("Free Evolution Times [ns]")
    ax.legend()
    if ax == None:
        return fig

def plot_multiT2( data, rep, time ):
    """
    data shape ( 2, N, M )
    2 is I,Q
    N is rep
    M is time
    """
    idata = data[0]
    qdata = data[1]
    zdata = idata +1j*qdata
    
    T2_list = []
    for i in range(len(rep)):
        T2 = T2_fitting(time,idata[i])
        T2_list.append(T2)
    T2_list = np.array(T2_list)
    print(f"T2 = {np.mean(T2_list)}")
    fig, ax = plt.subplots(2)
    ax[0].set_title('I signal')
    ax[0].pcolormesh( time, rep, idata, cmap='RdBu')# , vmin=z_min, vmax=z_max)
    ax[1].set_title('Q signal')
    ax[1].pcolormesh( time, rep, qdata, cmap='RdBu')# , vmin=z_min, vmax=z_max)
    return fig

