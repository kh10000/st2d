import os
import pickle

import numpy as np
from ase import io

from st2d.lib._libacsf import lib, ffi
from st2d.descriptors.basis import Basis
from st2d.descriptors.basis import _gen_2Darray_for_ffi

class ACSF(Basis):
    def __init__(self, input="./input.yaml"):
        self.basis_type = "acsf"
        super().__init__(input)
    

    def _read_params(self, filename):
        params_i = list()
        params_d = list()
        with open(filename, "r") as fil:
            for line in fil:
                tmp = line.split()
                params_i += [list(map(int,   tmp[:3]))]
                params_d += [list(map(float, tmp[3:]))]

        params_i = np.asarray(params_i, dtype=np.intc,    order="C")
        params_d = np.asarray(params_d, dtype=np.float64, order="C")

        return params_i, params_d


    def generate(self):
        """
        Generate descriptors from str_list.
        data*.pickle contains the descriptor vectors.
        """
        
        comm = super()._get_comm()
        
        if comm.rank == 0:
            print("Start the computation of descriptors.")
            pickle_list = open(self.inputs["pickle_list"], "w")

            super()._make_data_dir()
        
        # Get structure list to calculate
        structures = super()._parse_strlist(self.inputs["structure_list"])

        # Get parameter list for each atom types
        params_set = dict()
        for item in self.inputs["atom_types"]:    
            params_fil_name = self.inputs["params"] + "_" + item
            params_set[item] = dict()
            
            params_set[item]["i"], params_set[item]["d"] = \
                self._read_params(params_fil_name)
            
            params_set[item]["ip"] = _gen_2Darray_for_ffi(params_set[item]["i"], ffi, "int")
            params_set[item]["dp"] = _gen_2Darray_for_ffi(params_set[item]["d"], ffi)
            params_set[item]["total"] = np.concatenate((params_set[item]["i"], params_set[item]["d"]), axis=1)
            params_set[item]["num"] = len(params_set[item]["total"])

        data_idx = 1
        for structure in structures:
            atoms = io.read(structure)
            cell  = np.copy(atoms.cell, order="C")
            scale = np.copy(atoms.get_scaled_positions(), order="C")
            cart  = np.copy(np.dot(atoms.get_scaled_positions(), atoms.cell), order="C")

            symbols  = np.array(atoms.get_chemical_symbols())
            atom_num = len(symbols)
            atom_i   = np.zeros([len(symbols)], dtype=np.intc, order="C")
            type_num = dict()
            type_idx = dict()
            for j,jtem in enumerate(self.inputs["atom_types"]):
                tmp = symbols==jtem
                atom_i[tmp] = j+1
                type_num[jtem] = np.sum(tmp).astype(np.int64)
                type_idx[jtem] = np.arange(atom_num)[tmp]
            atom_i_p = ffi.cast("int *", atom_i.ctypes.data)

            res = dict()
            res["x"] = dict()

            cart_p  = _gen_2Darray_for_ffi(cart,  ffi)
            scale_p = _gen_2Darray_for_ffi(scale, ffi)
            cell_p  = _gen_2Darray_for_ffi(cell,  ffi)

            for j,jtem in enumerate(self.inputs["atom_types"]):
                q = type_num[jtem] // comm.size
                r = type_num[jtem] %  comm.size

                begin = comm.rank * q + min(comm.rank, r)
                end = begin + q
                if r > comm.rank:
                    end += 1
                
                cal_atoms = np.asarray(type_idx[jtem][begin:end], dtype=np.intc, order="C")
                cal_num = len(cal_atoms)
                cal_atoms_p = ffi.cast("int *", cal_atoms.ctypes.data)

                x = np.zeros([cal_num, params_set[jtem]["num"]], dtype=np.float64, order="C")
                x_p = _gen_2Darray_for_ffi(x, ffi)

                errno = lib.calculate_sf(cell_p, cart_p, scale_p, \
                                         atom_i_p, atom_num, cal_atoms_p, cal_num, \
                                         params_set[jtem]["ip"], params_set[jtem]["dp"], params_set[jtem]["num"], \
                                         x_p)
                
                comm.barrier()
                errnos = comm.gather(errno)
                errnos = comm.bcast(errnos)
                for ernno in errnos:
                    if errno == 1:
                        err = "Not implemented symmetry function type."
                        if comm.rank == 0:
                            print("\nError: {:}\n".format(err))
                        raise NotImplementedError(err)
                    elif ernno == 2:
                        err = "Zeta in G4/G5 must be greater or equal to 1.0."
                        if comm.rank == 0:
                            print("\nError: {:}\n".format(err))
                        raise ValueError(err)
                    else:
                        assert ernno == 0
                
                if type_num[jtem] != 0:
                    res["x"][jtem] = comm.gather(x, root=0)
                    # res["x"][jtem] = np.array(comm.gather(x, root=0), dtype=object)
                    if comm.rank == 0:
                        # res["x"][jtem] = np.concatenate(res["x"][jtem], axis=0).\
                        #             reshape([type_num[jtem], params_set[jtem]["num"]])
                        res["x"][jtem] = np.concatenate(res["x"][jtem], axis=0)
                        
                else:
                    res["x"][jtem] = np.zeros([0, params_set[jtem]["num"]])
                
            if comm.rank == 0:
                tmp_filename = os.path.join(self.data_dir, "data{}.pickle".format(data_idx))

                with open(tmp_filename, "wb") as fil:
                    pickle.dump(res, fil)
                    
                pickle_list.write("{}:{}\n".format(data_idx, tmp_filename))
                tmp_endfile = tmp_filename

                data_idx += 1
            
        if comm.rank == 0:
            print(": ~{}\n".format(tmp_endfile))
            print("Finish!")
            pickle_list.close()