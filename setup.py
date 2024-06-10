"""
Setup script for the UppASD Python wrapper.
"""
import os
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    """
    A custom extension class for CMake-based extensions.

    This class extends the `Extension` class and provides additional functionality
    for CMake-based extensions. It allows specifying the CMakeLists.txt directory
    and provides an abstraction for building CMake-based extensions.

    Args:
        name (str): The name of the extension.
        cmake_lists_dir (str, optional): The directory containing the CMakeLists.txt file.
            Defaults to the current directory.
        **kwa: Additional keyword arguments to be passed to the `Extension` class.

    Attributes:
        cmake_lists_dir (str): The absolute path to the CMakeLists.txt directory.

    """

    def __init__(self, name, cmake_lists_dir=".", **kwa):
        Extension.__init__(self, name, sources=[], **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)


class CmakeBuildExt(build_ext):
    """
    Custom build_ext class for building CMake extensions.

    This class extends the build_ext class from setuptools to provide
    custom build behavior for CMake extensions. It ensures that CMake
    is present and working, and then builds the extensions using CMake.

    Attributes:
        build_temp (str): The directory where the build files will be placed.
    """

    def build_extensions(self):

        # Ensure that CMake is present and working
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError("Cannot find CMake executable")

        for ext in self.extensions:
            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            cfg = (
                "Debug"
                if os.environ.get("DISPTOOLS_DEBUG", "OFF") == "ON"
                else "Release"
            )
            cmake_args = [
                f"-DCMAKE_BUILD_TYPE={cfg}",
                # Ask CMake to place the resulting library in the directory
                # containing the extension
                f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}",
                # Hint CMake to use the same Python executable that
                # is launching the build, prevents possible mismatching if
                # multiple versions of Python are installed
                # f"-DPython3_ROOT_DIR={sys.exec_prefix}",
                "-DPython3_FIND_STRATEGY=LOCATION",
                "-DBUILD_PYTHON=ON",
                "-DUSE_MKL=OFF",
                "-GNinja",
            ]

            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)

            # Config
            subprocess.check_call(
                ["cmake", ext.cmake_lists_dir] + cmake_args, cwd=self.build_temp
            )

            # Build
            subprocess.check_call(
                ["cmake", "--build", ".", "--parallel", "--config", cfg],
                cwd=self.build_temp,
            )

setup(
    # name = 'uppasd',
    # version = '1.0.1',
    # description = 'An UppASD Python wrapper',
    # url = 'https://github.com/UppASD/UppASD',
    # author = 'UppASD group',
    # author_email = 'uppasd@physics.uu.se',
    # license = 'GPLv2',
    cmdclass={"build_ext": CmakeBuildExt},
    # include_package_data = True,
    # packages=['uppasd'],
    # package_dir={'uppasd': 'uppasd'},
    ext_modules=[CMakeExtension(name="_uppasd")],
    # scripts=['scripts/uppasd','scripts/uppasd_interactive'],
    # install_requires=['numpy>1.19']
)
