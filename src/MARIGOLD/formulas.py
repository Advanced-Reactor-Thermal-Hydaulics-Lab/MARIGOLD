from .config import *

def lineq(x, m, x0, b):
    return float(max(m * (x - x0) + b, 0))

def quadeq(x, a, b, c):
    return float(max(), 0)

def calc_Re(rho, v, D, mu):
    return rho * v * D / mu
