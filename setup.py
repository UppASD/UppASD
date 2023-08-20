import os
import sys
import json
import platform
import shutil
import glob
from setuptools import setup, find_packages
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sysconfig



class CMakeExtension(Extension):
    #def __init__(self, name, cmake_lists_dir='../', **kwa):
    def __init__(self, name, cmake_lists_dir='.', **kwa):
        Extension.__init__(self, name, sources=[], **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)

class cmake_build_ext(build_ext):
    def build_extensions(self):

        import subprocess

        # Ensure that CMake is present and working
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError('Cannot find CMake executable')

        for ext in self.extensions:

            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            cfg = 'Debug' if os.environ.get('DISPTOOLS_DEBUG','OFF') == 'ON' else 'Release'
            python_inc_dir = sysconfig.get_path('include')
            python_lib_dir = sysconfig.get_config_var('LIBDIR')


            cmake_args = [
                '-DCMAKE_BUILD_TYPE=%s' % cfg,
                # Ask CMake to place the resulting library in the directory
                # containing the extension
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                # Other intermediate static libraries are placed in a
                # temporary build directory instead
                #'-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                #'-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), self.build_temp),
                # Hint CMake to use the same Python executable that
                # is launching the build, prevents possible mismatching if
                # multiple versions of Python are installed
                '-DPython3_ROOT_DIR={}'.format(sys.exec_prefix),
                '-DPython3_FIND_STRATEGY=LOCATION',
                #'-DCMAKE_Fortran_COMPILER=gfortran',
                '-DBUILD_PYTHON=ON',
                '-GNinja',
                #'-DPYTHON_INCLUDE_DIR={}'.format(python_inc_dir),
                #'-DPYTHON_LIBRARY={}'.format(python_lib_dir),
                #'-DMKL_INTERFACE_FULL=gf_lp64',
                #'-DMKL_THREADING=gnu_thread',
                #'-DLAPACK="-framework Accelerate"',
                #'-DBLAS="-framework Accelerate"',
            ]

            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)

            # Config
            subprocess.check_call(['cmake', ext.cmake_lists_dir] + cmake_args,
                                  cwd=self.build_temp)

            # Build
            subprocess.check_call(['cmake', '--build', '.','--parallel', '--config', cfg],
                                  cwd=self.build_temp)

            #src_file=glob.glob('./'+self.build_temp+'/_uppasd.*.*')
            #lib_path=self.build_temp.replace('temp','lib') #+'/uppasd/'
            #if not os.path.exists(lib_path):
            #    os.makedirs(lib_path)
            #
            #shutil.copy2(src_file[0],'uppasd/')
            #shutil.copy2(src_file[0],lib_path)
            


##### Environment flag needed for appending library flags to f2py
#os.environ['NPY_DISTUTILS_APPEND_FLAGS'] = "1"
##### Special flags to 
#if (platform.system()=='Darwin'):
#    os.environ['LDFLAGS'] = "-framework Accelerate"
#elif (platform.system()=='Linux'):
#    os.environ['LDFLAGS'] = "-fopenmp"

setup(
        name = 'uppasd',
        version = '1.0.1',
        description = 'An UppASD Python wrapper',
        url = 'https://github.com/UppASD/UppASD',
        author = 'UppASD group',
        author_email = 'uppasd@physics.uu.se',
        license = 'GPLv2',
        cmdclass={'build_ext':cmake_build_ext},
        include_package_data = True, 
        packages=['uppasd'],
        package_dir={'uppasd': 'uppasd'},
        ext_modules=[CMakeExtension(name='_uppasd')],
        scripts=['scripts/uppasd','scripts/uppasd_interactive'],
        install_requires=['numpy>1.19']
        )
