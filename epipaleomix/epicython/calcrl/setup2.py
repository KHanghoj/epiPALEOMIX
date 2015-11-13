from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
## import numpy as np                           # <---- New line
import pysam
## command: python setup2.py build_ext -i


from Cython.Build import cythonize
extension = [Extension("getminmax",
                       sources=["getminmax.pyx"],
                       include_dirs=pysam.get_include())
]

setup(
      name = 'getminmax',
      ext_modules = cythonize(extension)
)

# ext_modules = [Extension("getminmax",
#                          sources=["getminmax.pyx"],
#                          include_dirs=pysam.get_include())] #,
#                          # define_macros=pysam.get_defines(),
#                          # extra_link_args=pysam.get_libraries())]


# setup(
#       name = 'getminmax',
#       cmdclass = {'build_ext': build_ext},
#       ext_modules = ext_modules
# )


### if need to add several modules to the same nucleopy.pyx file do like this:
# ext_modules = [Extension("nucleoPynew",
#                          sources=["getminmax.pyx"],
#                          include_dirs=[np.get_include()]),
#                Extension("nucleoPynew",
#                          sources=["getminmax.pyx"],
#                          include_dirs=[pysam.get_include()],
#                          define_macros=pysam.get_defines(),
#                          extra_link_args=pysam.get_libraries())
# ]

# setup(
#       name = 'nucleoPynew',
#       cmdclass = {'build_ext': build_ext},
#       ext_modules = ext_modules
# )


## http://docs.cython.org/src/reference/compilation.html
# from distutils.core import setup
# from distutils.extension import Extension
# from Cython.Build import cythonize

# extensions = [
#     Extension("primes", ["primes.pyx"],
#               include_dirs = [...],
#               libraries = [...],
#               library_dirs = [...]),
#     # Everything but primes.pyx is included here.
#     Extension("*", ["*.pyx"],
#               include_dirs = [...],
#               libraries = [...],
#               library_dirs = [...]),
# ]
# setup(
#     name = "My hello app",
#     ext_modules = cythonize(extensions),
# )
