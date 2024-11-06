from qm.qua import *
from qm.QuantumMachinesManager import QuantumMachinesManager
from qm import SimulationConfig
from qualang_tools.results import progress_counter, fetching_tool
from qualang_tools.plot import interrupt_on_close
from qualang_tools.loops import from_array
from exp.RO_macros import multiRO_declare, multiRO_measurement, multiRO_pre_save
import matplotlib.pyplot as plt
import warnings
# from common_fitting_func import *
warnings.filterwarnings("ignore")
import exp.config_par as gc
from qualang_tools.units import unit
u = unit(coerce_to_integer=True)

import xarray as xr
import time
from exp.QMMeasurement import QMMeasurement

class DetunedRabiFluxPulse( QMMeasurement ):
    """
    Parameters:\n
    freq_span:\n
        a tuple (upper, lower) Unit in MHz, \n
    amp_max_ratio: \n

    amp_scale: \n
        lin or log \n
    Return: xarray dataset
        coords : frequency, amp_ratio
    """
    def __init__( self, config, qmm: QuantumMachinesManager ):
        super().__init__( config, qmm )

        self.ro_elements = ["q4_ro"]
        self.z_elements = ["q4_z"]
        self.xy_elements = ["q4_xy"]

        self.initializer = None
        

        self.duration = 400
        self.time_resolution = 80
        self.pad_zeros = ( 80, 0 )

        self.amp_modify = 0.2

        self.freq_range = ( -200, 50 )
        self.freq_resolution = 100

        

    def _get_qua_program( self ):
        
        self.qua_cc_pi_timing = self._lin_cc_array()
        # print(self.qua_cc_pi_timing)
        self.qua_freqs = self._lin_freq_array()
        self.qua_cc_duration = self.duration//4

        self._attribute_config()
        
        with program() as multi_res_spec_vs_amp:
        
            iqdata_stream = multiRO_declare( self.ro_elements )
            n = declare(int)
            n_st = declare_stream()
            df = declare(int)
            cc = declare(int)

            with for_(n, 0, n < self.shot_num, n + 1):
                with for_(*from_array(cc, self.qua_cc_pi_timing)):
                    
                    with for_(*from_array(df, self.qua_freqs)):
                        # Initialization
                        if self.initializer is None:
                            wait(1*u.us, self.ro_elements)
                        else:
                            try:
                                self.initializer[0](*self.initializer[1])
                            except:
                                print("initializer didn't work!")
                                wait(1*u.us, self.ro_elements)

                        # Operation
                        align()
                        for z in self.z_elements:
                            wait( (16+self.pad_zeros[0])//4,z)
                            play("const"*amp(self.amp_modify), z, duration=self.qua_cc_duration)
                            play("const"*amp(-0.04), "q8_z", duration=self.qua_cc_duration)
                        for i, xy in enumerate(self.xy_elements):
                            update_frequency( xy, self.ref_xy_IF[i] +df )

                            wait( cc, xy )
                            play("x180", xy )
                            play("x180", "q3_xy")
                        align()
                        # Readout
                        multiRO_measurement( iqdata_stream, self.ro_elements, weights='rotated_' )

                save(n, n_st)

            with stream_processing():
                n_st.save("iteration")
                # Cast the data into a 2D matrix, average the 2D matrices together and store the results on the OPX processor
                # NOTE that the buffering goes from the most inner loop (left) to the most outer one (right)
                multiRO_pre_save( iqdata_stream, self.ro_elements, (len(self.qua_cc_pi_timing), len(self.qua_freqs)))



        return multi_res_spec_vs_amp
    
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
            output_data[r_name] = ( ["mixer","pi_time","frequency"],
                                np.array([ self.fetch_data[r_idx*2], self.fetch_data[r_idx*2+1]]) )

        freqs_mhz = self.qua_freqs/1e6
        pi_timing = self.qua_cc_pi_timing*4

        dataset = xr.Dataset(
            output_data,
            coords={ "mixer":np.array(["I","Q"]), "frequency": freqs_mhz, "pi_timing": pi_timing }
        )

        self._attribute_config()
        dataset.attrs["ro_LO"] = self.ref_ro_LO
        dataset.attrs["ro_IF"] = self.ref_ro_IF
        dataset.attrs["xy_LO"] = self.ref_xy_LO
        dataset.attrs["xy_IF"] = self.ref_xy_IF
        dataset.attrs["z_offset"] = self.z_offset

        dataset.attrs["z_amp"] = self.z_amp

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
        self.z_amp = []
        for z in self.z_elements:
            self.z_offset.append( gc.get_offset(z, self.config ))
            self.z_amp.append(gc.get_const_wf(z, self.config ))

    def _lin_freq_array( self ):

        freq_r1_qua = self.freq_range[0] * u.MHz
        freq_r2_qua = self.freq_range[1] * u.MHz
        freq_resolution_qua = self.freq_resolution * u.MHz
        freqs_qua = np.arange(freq_r1_qua,freq_r2_qua,freq_resolution_qua )
        
        return freqs_qua

    def _lin_cc_array( self ):
         
        qua_cc_resolution = self.time_resolution //4
        if qua_cc_resolution < 1: qua_cc_resolution = 1

        cc_array = np.arange( 16//4, (16+self.duration+self.pad_zeros[0]+self.pad_zeros[1])//4, qua_cc_resolution)
        return cc_array
    

