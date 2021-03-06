=========
Chebyshev
=========
In this section, the procedure for generating Chebyshev descriptors[1] is explained.
The minimum input files required for the calculation are ``input.yaml``, ``str_list``, and hyper-parameter setting files.
The example of ``input.yaml`` in Ti-O-Sr system is shown below.

.. code-block:: text

    structure_list: ./str_list
    data_dir: ./data
    pickle_list: ./pickle_list
    atom_types:
     - Sr
     - Ti
     - O
    atom_weights:
     - -1
     - 0
     - 1
    params: params

Here, ``structure_list`` specifies the file containing the paths to the structure files to be calculated. The example is shown bellow. 
Structure files can be in any format as long as they can be recognized by **ase.io.read()**.

.. code-block:: text

    ../structure_data/structure_file_0
    ../structure_data/structure_file_1
    ../structure_data/structure_file_2

``data_dir`` specifies the path to the directory to which output files (pickle files) are stored, and ``pickle_list`` specifies the file where the name of output pickle files are recorded.
All the atomic species and the corresponding weights must be listed in ``atom_types`` and ``atom_weights``.

Hyper-parameters for each atom type is parsed acccording to the value of ``params``. Specifically, file with the name of "{value}_{atom type}" is parsed.
For example, if "params" is set in ``params``, three files, ``params_Sr``, ``params_Ti``, and ``params_O`` should be prepared. An example of ``params_X`` is as follows.

.. code-block:: text

    5 7.0 # RDF
    5 7.0 # weighted-RDF
    5 7.0 # ADF
    5 7.0 # weighted-ADF

The values in the first column are the expansion order, and that in the second column are cutoff radius.
The first through fourth lines correspond to radial distribution function (RDF), weighted RDF, angular distribution function (ADF), and weighted ADF, respectively.

Once ``input.yaml``, ``str_list``, and ``params_X`` are prepared, the descriptors can be generated by runnig the following script.

.. code-block:: python

    from st2d.descriptors.chebyshev import Chebyshev

    chebyshev = Chebyshev("input.yaml")
    chebyshev.generate()

Reference
---------
| [1] N. Artrith *et al.*, Phys. Rev. B **96**, 014112 (2017).