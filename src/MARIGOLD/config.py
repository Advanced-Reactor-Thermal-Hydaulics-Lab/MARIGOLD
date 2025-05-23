import os
import numpy as np
from datetime import datetime
import pickle
from scipy import integrate
from scipy import interpolate
from copy import deepcopy
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib as mpl
import openpyxl as op
import warnings
import re

debug = False
debugFID = None

__all__ = ['os', 'np', 'datetime', 'pickle', 'integrate', 'deepcopy', 'cm', 'plt', 'tri', 'op', 're', 'debug', 'warnings','debugFID']