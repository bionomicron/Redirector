#!/usr/bin/pythong
'''
Setup script for Redirector
@author: Graham Rockwell
20130202
'''
import os
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

def read(fname):
    '''
    Utility function to read the files file.
    '''
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


#Description = open('README').read()

#License = open('LICENSE').read()

#Version = open('VERSION').read().strip()

setup(
    name = "Redirector",
    version = "0.1",
    author = "Graham Rockwell",
    author_email = "grockwell@receptor.med.harvard.edu",
    description = ("Framework for designing metabolic engineering targets for optimization of cellular factories"),
    #long_description=Description,
    #license = License,
    keywords = ["Optimization", "Linear Programming", "Flux Balance Analysis", "Metabolic Objective"],
    url = "http://packages.python.org/an_example_pypi_project",
    classifiers=[
        "Development Status :: Alpha",
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Metabolic Modeling',        
    ],
    packages=find_packages(),
    install_requires=['PuLP==1.5.0','gmpy>=1.10'],
    entry_points = ("""
    [console_scripts]
    redirector = core.analysis.Redirector:main_function
    """
    ),
)