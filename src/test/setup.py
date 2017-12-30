from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("C_Compare", ["C_Compare.pyx"],include_dirs=[numpy.get_include()])]

setup(
  name = 'Compare',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)




#python setup.py build_ext --inplace