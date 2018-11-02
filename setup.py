from distutils.core import setup
from Cython.Build import cythonize

extra_compile_args=["-O3"]

setup(
    #ext_modules = cythonize('global_coupling.pyx')
    ext_modules = cythonize('global_coupling_largeN.pyx')
    #ext_modules = cythonize('global_coupling_largeN_try.pyx')
)
