# -*- coding: utf-8 -*-
from setuptools import setup, find_packages
import sys

setup(
    name="pepFoot",
    version="1.2.1",
    author="Jedd Bellamy-Carter",
    description="Graphical user interface application for semi-automated processing of protein footprinting data",
    license="LGPL-3.0",
    url="https://github.com/jbellamycarter/pepfoot",
    packages=find_packages(),
    install_requires=['PyQt5>=5.11',
                      'PyQtWebEngine',
                      'matplotlib',
                      'numpy',
                      'scipy',
                      'pyteomics',
                      'h5py',
                      'requests'],
    package_data={'pepFoot':['ngl.js','gui/*.png','COPYING*', '*.md']},
    entry_points={'gui_scripts':['pepfoot = pepFoot.pepFootGui:run_app']},
    classifiers=["License :: OSI Approved :: LGPL-3.0 License"]
)
