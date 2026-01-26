import os
import sys
import platform
import shutil
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sysconfig


class CMakeExtension(Extension):
    # def __init__(self, name, cmake_lists_dir='../', **kwa):
    def __init__(self, name, cmake_lists_dir=".", **kwa):
        Extension.__init__(self, name, sources=[], **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)


class cmake_build_ext(build_ext):
    def build_extensions(self):
        import subprocess

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
            python_inc_dir = sysconfig.get_path("include")
            python_lib_dir = sysconfig.get_config_var("LIBDIR")

            cmake_args = [
                "-DCMAKE_BUILD_TYPE=%s" % cfg,
                # Place resulting library in the extension directory
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir),
                # Use the current Python
                "-DPython3_ROOT_DIR={}".format(sys.exec_prefix),
                "-DPython3_FIND_STRATEGY=LOCATION",
                "-DBUILD_PYTHON=ON",
                "-DUSE_MKL=OFF",
                # Ensure gfortran is used rather than f95 alias when available
                "-DCMAKE_Fortran_COMPILER=gfortran",
            ]

            # Select a generator robustly: prefer Ninja if available, otherwise let CMake default
            cmake_generator = os.environ.get("CMAKE_GENERATOR")
            if cmake_generator:
                cmake_args.append(f"-G{cmake_generator}")
            else:
                if shutil.which("ninja"):
                    cmake_args.append("-GNinja")
                # else: rely on CMake default (Unix Makefiles on Linux)

            # Only set macOS deployment target on macOS
            if platform.system() == "Darwin":
                # Default to macOS 12 if not specified via env
                mac_deploy_target = os.environ.get("CMAKE_OSX_DEPLOYMENT_TARGET", "12")
                cmake_args.append(f"-DCMAKE_OSX_DEPLOYMENT_TARGET={mac_deploy_target}")

            # BLAS/LAPACK vendor selection: allow override via env; don't force a default
            bla_vendor = os.environ.get("BLA_VENDOR")
            if bla_vendor:
                cmake_args.append(f"-DBLA_VENDOR={bla_vendor}")
            # Otherwise, rely on CMake's FindBLAS/FindLAPACK to resolve the
            # available implementation (conda libblas on Binder/Colab, Accelerate on macOS, etc.).

            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)

            # Use an out-of-tree build: explicitly set source and build directories
            # This prevents CMake from doing an in-tree configure in the source tree.
            subprocess.check_call(
                ["cmake", "-S", ext.cmake_lists_dir, "-B", self.build_temp] + cmake_args
            )

            # Build using the build directory
            subprocess.check_call(
                ["cmake", "--build", self.build_temp, "--parallel", "--config", cfg]
            )


setup(
    # name="uppasd",
    # version="1.0.1",
    # description="An UppASD Python wrapper",
    # url="https://github.com/UppASD/UppASD",
    # author="UppASD group",
    # author_email="uppasd@physics.uu.se",
    # license="GPLv2",
    cmdclass={"build_ext": cmake_build_ext},
    # include_package_data=True,
    # packages=["uppasd"],
    # package_dir={"uppasd": "uppasd"},
    ext_modules=[CMakeExtension(name="_uppasd")],
    # scripts=["scripts/uppasd", "scripts/uppasd_interactive"],
    # install_requires=["numpy>1.19"],
)
