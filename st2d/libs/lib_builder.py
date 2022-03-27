import cffi

# acsf
ffibuilder_sf = cffi.FFI()
ffibuilder_sf.cdef(
    """int calculate_sf(double **, double **, double **,
                                    int *, int, int*, int,
                                    int**, double **, int,
                                    double**);"""
)
ffibuilder_sf.set_source(
    "st2d.lib._libacsf",
    '#include "calculate_sf.h"',
    sources=[
        "st2d/lib/calculate_sf.cpp",
    ],
    source_extension=".cpp",
    include_dirs=["st2d/lib/"],)

# chebyshev
ffibuilder_ch = cffi.FFI()
ffibuilder_ch.cdef(
    """int calculate_chebyshev(double **, double **, double **,
                                    int *, double *, int, int*, int,
                                    int*, double *, int,
                                    double**);"""
)
ffibuilder_ch.set_source(
    "st2d.lib._libchebyshev",
    '#include "calculate_chebyshev.h"',
    sources=[
        "st2d/lib/calculate_chebyshev.cpp",
    ],
    source_extension=".cpp",
    include_dirs=["st2d/lib/"],)



if __name__ == "__main__":
    ffibuilder_sf.compile(verbose=True)
    ffibuilder_ch.compile(verbose=True)
