from setuptools import setup, find_packages
import json

setup(
        name = 'asd_gui',
        version = '0.9',
        description = 'The UppASD GUI',
        url = 'https://github.com/UppASD/UppASD',
        author = 'Jonathan P. Chico',
        author_email = 'uppasd@physics.uu.se',
        license = 'GPLv2',
        packages=find_packages(),
        package_dir={'ASD_GUI.UI': 'ASD_GUI/UI'},
        package_data={'ASD_GUI': ['*.ui']},
        include_package_data = True,
        scripts=['bin/asd_gui'],
        install_requires=[
            'matplotlib',
            'numpy>=1.19.0',
            'pandas',
            'PyQt6',
            'PyYAML',
            'setuptools',
            'vtk>=9.0.0',
            'scipy',
            ]
        )

