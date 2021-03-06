import os, sys

import yaml

from st2d.utils.mpiclass import MPI4PY, DummyMPI

class Descriptor():
    def __init__(self, input="./input.yaml"):
        self.comm = None

        self.input_fil = input
        self.default_inputs = {
            "structure_list":"./str_list",
            "data_dir":"./data",
            "pickle_list":"./pickle_list",
            "atom_types":[],
            "params":"params"
        }
        self.inputs = self.default_inputs
        if self.descriptor_type == "Chebyshev":
            self.inputs["atom_weights"] = []
        self._set_inputs(self.input_fil)


    def _set_inputs(self, infile):
        with open(infile, "r") as yml:
            input = yaml.safe_load(yml)
        for key in self.inputs:
            try:
                self.inputs[key] = input[key]
            except:
                print(f"Warning: {key} was not set.")
                pass
    

    def _make_data_dir(self):
        self.data_dir = self.inputs["data_dir"]
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
    

    def _parse_strlist(self, str_list):
        structures = []
        with open(str_list, "r") as fil:
            for line in fil:
                line = line.strip()
                structures.append(line)
        
        return structures
    
    
    def _get_comm(self):
        if self.comm is None and "mpi4py" not in sys.modules:
            try:
                import mpi4py
            except ImportError:
                self.comm = DummyMPI()
            else:
                self.comm = MPI4PY()
                
        elif self.comm is None:
            self.comm = MPI4PY()
        
        else:
            self.comm = self.comm

        return self.comm


def _gen_2Darray_for_ffi(arr, ffi, cdata="double"):
    """
    Function to generate 2D pointer for cffi  
    """
    shape = arr.shape
    arr_p = ffi.new(cdata + " *[%d]" % shape[0])
    for i in range(shape[0]):
        arr_p[i] = ffi.cast(cdata + " *", arr[i].ctypes.data)

    return arr_p    
