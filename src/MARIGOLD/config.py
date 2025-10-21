import inspect
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import openpyxl as op
import os
import pickle
import re
import subprocess
import warnings
import xlrd

from copy import copy
from copy import deepcopy
from datetime import datetime
from matplotlib import cm
from scipy import integrate
from scipy import interpolate
from scipy.optimize import minimize
from shutil import copy2
from subprocess import run

debug = False
debugFID = None

# Automatically generate __all__ from current namespace
__all__ = [
    name for name in dir()
    if not name.startswith('_')  # skip private names
]
