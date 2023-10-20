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
from openpyxl import load_workbook

debug = False
debugFID = None

__all__ = ['os', 'np', 'datetime', 'pickle', 'integrate', 'deepcopy', 'cm', 'plt', 'tri', 'load_workbook', 'debug', 'debugFID']