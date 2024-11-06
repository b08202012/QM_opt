
from qm.qua import *
from qm.QuantumMachinesManager import QuantumMachinesManager
from qualang_tools.units import unit
u = unit(coerce_to_integer=True)
from qualang_tools.loops import from_array
import warnings

from exp.RO_macros import multiRO_declare, multiRO_measurement, multiRO_pre_save

warnings.filterwarnings("ignore")
import exp.config_par as gc
import xarray as xr
from exp.QMMeasurement import QMMeasurement

import numpy as np

class ROFreqSweep( QMMeasurement ):
    """
    Parameters:
    Search cavities with the given IF range along the given ro_element's LO, (LO+freq_range[0],LO+freq_range[1]) .\n

    freq_range:\n
    is a tuple ( upper, lower ), Unit in MHz.\n
    resolution:\n
    unit in MHz.\n
    ro_element: ["q1_ro"], temporarily support only 1 element in the list.\n
    initializer: from `initializer(paras,mode='depletion')`, and use paras return from `Circuit_info.give_depletion_time_for()`  
    Return: \n

    """
    def __init__( self, config, qmm: QuantumMachinesManager ):
        super().__init__( config, qmm )

        self.ro_elements = ["q0_ro"]
        self.initializer = None
        
        self.freq_range = ( -300, 300 )
        self.resolution = 1.

    

    def _get_qua_program( self ):

        self.frequencies_qua = self._lin_freq_array( )
        self._attribute_config()

        with program() as resonator_spec:

            f = declare(int)  # QUA variable for the readout frequency --> Hz int 32 up to 2^32
            iqdata_stream = multiRO_declare( self.ro_elements )
            n = declare(int)
            n_st = declare_stream()

            with for_(n, 0, n < self.shot_num, n + 1):  # QUA for_ loop for averaging
                with for_(*from_array(f, self.frequencies_qua)):  # QUA for_ loop for sweeping the frequency
                    # Initialization
                    # Wait for the resonator to deplete
                    if self.initializer is None:
                        wait(10 * u.us, self.ro_elements[0])
                    else:
                        try:
                            self.initializer[0](*self.initializer[1])
                        except:
                            print("initializer didn't work!")
                            wait(1 * u.us, self.ro_elements[0]) 
                    # Operation
                    update_frequency( self.ro_elements[0], f)
                    # Readout
                    multiRO_measurement( iqdata_stream, self.ro_elements, weights="rotated_") 
                # Save the averaging iteration to get the progress bar
                save(n, n_st)

            with stream_processing():
                # Cast the data into a 1D vector, average the 1D vectors together and store the results on the OPX processor
                multiRO_pre_save( iqdata_stream, self.ro_elements, (len(self.frequencies_qua),))
                n_st.save("iteration")

        return resonator_spec

    def _get_fetch_data_list( self ):
        return [f"{self.ro_elements[0]}_I", f"{self.ro_elements[0]}_Q", "iteration"]

    def _data_formation( self ):

        frequencies_mhz = self.frequencies_qua/1e6 #  Unit in MHz
        output_data = np.array([self.fetch_data[0],self.fetch_data[1]])
        dataset = xr.Dataset(
            {
                self.ro_elements[0]: (["mixer","frequency"], output_data),
            },
            coords={"frequency": frequencies_mhz, "mixer":np.array(["I","Q"]) }
        )
        self._attribute_config()
        dataset.attrs["ro_LO"] = self.ref_ro_LO
        dataset.attrs["ro_IF"] = self.ref_ro_IF

        return dataset

    def _lin_freq_array( self ):
        return np.arange( self.freq_range[0]*u.MHz, self.freq_range[1]*u.MHz, self.resolution* u.MHz )
    

    def _attribute_config( self ):
        self.ref_ro_IF = []
        self.ref_ro_LO = []
        for r in self.ro_elements:
            self.ref_ro_IF.append(gc.get_IF(r, self.config))
            self.ref_ro_LO.append(gc.get_LO(r, self.config))


