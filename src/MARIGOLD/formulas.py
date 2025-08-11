from .config import *
import warnings

def calc_Re():
    
    return

def calc_Reb(cond):
    """TODO, not impelmented
    
    """
    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            Reb = 0
            midas_dict.update({'Reb': Reb})

