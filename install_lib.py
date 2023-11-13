"""
install_lib.py
Script install relevant libraries for USC workshop
Created on Tue Jun 20 16:42:31 2023
@author:Carlos Alberto Duran Villalobos-The University of Manchester
"""



import sys
import subprocess
#import conda.cli.python_api as Conda


def install_lib():   
    try:
        import numpy
    except ModuleNotFoundError:
        print("numpy not found")
        subprocess.check_call([sys.executable, '-m', 'conda', 'install', 'numpy'])
    try:
        import sklearn
    except ModuleNotFoundError:
        print("sklearn not found")
        subprocess.check_call([sys.executable, '-m', 'conda', 'install', 'sklearn'])
    try:
        import scipy
    except ModuleNotFoundError:
        print("scipy not found")
        subprocess.check_call([sys.executable, '-m', 'conda', 'install', 'scipy']) 
    try:
        import seaborn
    except ModuleNotFoundError:
        print("sb not found")
        subprocess.check_call([sys.executable, '-m', 'conda', 'install', 'seaborn'])
    try:
        import ipywidgets
    except ModuleNotFoundError:
        print("ipywidgets not found")
        subprocess.check_call([sys.executable, '-m', 'conda', 'install', 'ipywidgets']) 
    try:
        import matplotlib
    except ModuleNotFoundError:
        print("plt not found")
        subprocess.check_call([sys.executable, '-m', 'conda', 'install', 'matplotlib']) 
    try:
        import widgetsnbextension
    except ModuleNotFoundError:
        print("nb not found")
        subprocess.check_call([sys.executable, '-m', 'conda', 'install', 'widgetsnbextension']) 
