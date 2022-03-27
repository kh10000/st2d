from setuptools import setup

setup(
    cffi_modules = [
        "st2d/lib/lib_builder.py:ffibuilder_sf",
        "st2d/lib/lib_builder.py:ffibuilder_ch",
    ],
)