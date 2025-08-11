from .config import *
from .plot_utils import *
from scipy import interpolate
from scipy.optimize import minimize
import warnings

def approx_vf(cond, n=7, overwrite_vf = False) -> None:
    """Method for approximating :math:`v_{f}` with power-law relation. 

    .. math:: v_{f, approx} = \\frac{(n+1)(2n+1)}{ (2n^{2})}  (j_{f} / (1- \\langle \\alpha \\rangle))  (1 - |r^{*}|)^{1/n}
    
    **Args**:
    
        - ``n``: power. Defaults to 7.
        - ``overwrite_vf``: will update ``'vf'`` even if it already exists. Defaults to False.
    
    **Returns**:
    
        - Area-averaged ``'vf_approx'``
        - Stores: 
            - ``'vf_approx'`` in ``midas_dict``
            - ``'vf'``, if it does not alread exist in ``midas_dict``
    """

    cond.mirror()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            vf_approx = (n+1)*(2*n+1) / (2*n*n) * (cond.jf / (1-cond.area_avg('alpha'))) * (1 - abs(rstar))**(1/n)
            if 'vf' in midas_dict.keys():
                if debug: print(f"approx_vf: data found for {angle}\t{rstar}", file=debugFID)
                if overwrite_vf:
                    midas_dict.update({'vf': vf_approx})
            else:
                midas_dict.update({'vf': vf_approx})
            
            midas_dict.update({'vf_approx': vf_approx})

    return cond.area_avg('vf_approx')

def approx_vg(cond, method = 'vr', n=7, update_ug1 = False) -> None:
    """Method for approximating :math:`v_{g}` with power-law relation. I don't think this makes sense
    
    **Args**:
    
        - ``method``: method to use for approximating :math:`v_{g}`. Defaults to 'vr'.
            - ``'power-law'``, pretty stupid for bubbly flows, not recommended

            .. math:: v_{g, approx} = \\frac{(n+1)(2n+1)}{ (2n^{2})}  (j_{g} / \\langle \\alpha \\rangle)  (1 - |r^{*}|)^{1/n}

            - ``'vrmodel'``

            .. math:: v_{g, approx} = v_{f} + v_{r, model}

            - ``'vr'``

            .. math:: v_{g, approx} = v_{f} + v_{r}

        - ``n``: power. Defaults to 7.
        - ``update_ug1``: overwrite ``ug1`` in ``midas_dict``. Defaults to False.
    
    **Raises**:
    
        - ``ValueError``: If unknown method selected
    
    **Returns**:
    
        - Area-averaged ``'vg_approx'``
        - Stores:
            - ``'vg_approx'`` in ``midas_dict``
            - ``'ug1'``, if it does not alread exist in midas_dict
    """

    cond.mirror()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            if method == 'power-law':
                vg_approx = (n+1)*(2*n+1) / (2*n*n) * (cond.jgloc / (cond.area_avg('alpha'))) * (1 - abs(rstar))**(1/n)
            elif method == 'vrmodel':
                if 'vr_model' not in midas_dict.keys():
                    cond.calc_vr_model()
                vg_approx = midas_dict['vf_approx'] + midas_dict['vr_model']
            elif method == 'vr':
                vg_approx = midas_dict['vf'] + midas_dict['vr']
            elif type(method) is float:
                vg_approx = method
            else:
                raise ValueError(f"Unknown method for approximating vg: {method}")
            
            if update_ug1 and 'ug1' not in midas_dict.keys():
                midas_dict.update({'ug1': vg_approx})
            
            midas_dict.update({'vg_approx': vg_approx})

    return cond.area_avg('vg_approx')

def approx_vf_Kong(cond, n=7) -> None:
    """Not currently implemented, right now a 1/nth power law thing
    
    **Args**:
    
        - ``n``: power. Defaults to 7.
    """

    cond.mirror()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            vf_approx = (n+1)*(2*n+1) / (2*n*n) * (cond.jf / (1-cond.area_avg('alpha'))) * (1 - abs(rstar))**(1/n)
            midas_dict.update({'vf': vf_approx})

    return

def calc_vf_lee(cond, K=1):
    """Calculate :math:`v_{f}`, :math:`j_{f}`, :math:`v_{r}` based on Lee et al. (2002) equation

    Really a model proposed by Bosio and Malnes (1968)

    .. math:: v_{f} = \\frac{ 1 }{ \\sqrt{1 - \\alpha^{2} / 2} } * \\sqrt{ \\frac{ 2 \\Delta p }{K \\rho_{f}} }
    
    **Args**:
    
        - ``K``: momentum coefficient. Defaults to 1.
    
    **Raises**:
    
        - ``NotImplementedError``: If ``'delta_p'`` not in ``midas_dict``
    
    **Returns**:
    
        - Area-average vf_lee
        - Stores:
            - ``'vf_lee'`` in ``'midas_dict'``
            - ``'jf_lee'`` in ``'midas_dict'``
            - ``'vr_lee'`` in ``'midas_dict'``

    """

    cond.mirror()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            try:
                dp = midas_dict['delta_p'] * 6894.757
            except:
                raise NotImplementedError("Δp needed for cacluclation of vf_lee")
            
            vf_lee = 1 / np.sqrt(1 - midas_dict['alpha']**2/2) * np.sqrt( 2 * dp / (K * cond.rho_f))
            
            midas_dict.update({'vf_lee': vf_lee})
            midas_dict.update({'jf_lee': (1-midas_dict['alpha'])* vf_lee})
            
            vg = midas_dict['ug1']
            if vg == 0:
                vr_lee = 0
            else:
                vr_lee = vg - vf_lee
            midas_dict.update({'vr_lee':  vr_lee})

    return cond.area_avg('vf_lee')

def calc_vf_naive(cond):
    """Calculate :math:`v_{f}`, :math:`j_{f}`, :math:`v_{r}` based on single-phase Pitot-tube equation

    .. math:: v_{f, naive} = \\sqrt{ \\frac{2 \\Delta p}{\\rho_{f}} }        
    
    **Raises**:
    
        - ``NotImplementedError``: If ``'delta_p'`` not in ``midas_dict``
    
    **Returns**:
    
        - Area-average ``'vf_naive'``
        - Stores:
            - ``'vf_naive'`` in ``'midas_dict'``
            - ``'jf_naive'`` in ``'midas_dict'``
            - ``'vr_naive'`` in ``'midas_dict'``
    """

    cond.mirror()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            try:
                vf_naive = midas_dict['vf_naive']
            except KeyError:
                try:
                    dp = midas_dict['delta_p'] * 6894.757
                except:
                    raise NotImplementedError("Δp needed for cacluclation of vf_naive")
                
                vf_naive = np.sqrt( 2*dp / cond.rho_f)
            
            midas_dict.update({'vf_naive': vf_naive})
            midas_dict.update({'jf_naive': (1-midas_dict['alpha'])* vf_naive})
            
            vg = midas_dict['ug1']
            if vg == 0:
                vr_naive = 0
            else:
                vr_naive = vg - vf_naive
            midas_dict.update({'vr_naive':  vr_naive})

    return cond.area_avg('vf_naive')

def calc_vr(cond, method = None, quiet = False) -> None:
    """Calculate relative velocity

    .. math:: v_{r} = v_{g,1} - v_{f}

    Note that if vg = 0, then this method says vr = 0. This will happen when no data is present, such as in the bottom of the pipe in horizontal, when this is not necessarily true.
    
    **Args**:
    
        - ``method``: method to calculate :math:`v_r`. Defaults to None.

            - ``'approx'``, uses ``'vf_approx'`` in ``midas_dict``
            - ``'None'``, uses ``'vf'`` in ``midas_dict``

        - ``quiet``: warn if approximating. Defaults to False.
    
    **Returns**:
    
        - Area-average vr
        - Stores ``'vr'`` in ``midas_dict``
    """

    cond.mirror()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            try:
                if method == None:
                    vf = midas_dict['vf']
                elif method == 'approx':
                    cond.approx_vf()
                    vf = midas_dict['vf_approx']
            except:
                if not quiet:
                    print("Warning: Approximating vf in calculating vr, since no data found")
                    warn_approx = False
                cond.approx_vf()
                vf = midas_dict['vf']
            vg = midas_dict['ug1']

            if vg == 0: # should be the same as α = 0, could maybe switch this to that
                vr = 0 # this is an assumption, similar to void weighting
                
            else:
                vr = vg - vf

            try:
                if abs( midas_dict['vr'] - vr ) < 0.00001:
                    pass
                elif not quiet:
                    print(f"Warning: vr already present for {rstar}, {angle}°, but doesn't match subtraction. Will update and overwrite")
                    midas_dict.update({'vr': vr})
                elif quiet:    
                    midas_dict.update({'vr': vr})
            except:
                midas_dict.update({'vr': vr})
                

    return cond.area_avg('vr')

def calc_vr2(cond, warn_approx = True) -> None:
    """Method for calculating relative velocity based on ``'ug2'``

    .. math:: v_{r,2} = v_{g, 2} - v_{f}
    
    Note that if vg2 = 0, then this method says vr2 = 0. This will happen when no data is present,
    such as in the bottom of the pipe in horizontal, when this is not necessarily true

    **Args**:
    
        - ``warn_approx``: flag for warning if :math`v_r` is not found and function is approximating. Defaults to True.
    
    **Returns**:
    
        - Area-average vr2
        - Stores ``'vr2'`` in ``midas_dict``
    """

    cond.mirror()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            try:
                vf = midas_dict['vf']
            except:
                if warn_approx:
                    print("Warning: Approximating vf in calculating vr, since no data found")
                    warn_approx = False
                cond.approx_vf()
                vf = midas_dict['vf']
            vg = midas_dict['ug2']

            if vg == 0: # should be the same as α = 0, could maybe switch this to that
                vr = 0 # this is an assumption, similar to void weighting
                
            else:
                vr = vg - vf

            try:
                if abs( midas_dict['vr2'] - vr ) < 0.00001:
                    pass
                else:
                    print(f"Warning: vr2 already present for {rstar}, {angle}°, but doesn't match subtraction. Will update and overwrite")
                    midas_dict.update({'vr2': vr})
            except:
                midas_dict.update({'vr2': vr})
                

    return cond.area_avg('vr2')

def calc_vgj(cond, warn_approx = True) -> None:
    """Method for calculating :math:`V_{gj}`
    
    **Args**:
    
        - ``warn_approx``: flag for warning if :math:`V_{gj}` is not found and function is approximating. Defaults to True.
    
    **Returns**:
    
        - Void-weighted area-averaged Vgj
        - Stores:
            - ``'vgj'`` in ``midas_dict``
            - ``'j'`` in ``midas_dict``
            - ``'alpha_j'`` in ``midas_dict``
    """

    cond.mirror()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            try:
                dummy = midas_dict['vf']
            except:
                if warn_approx:
                    print("Warning: Approximating vf in calculating local j, since no data found")
                    warn_approx = False
                cond.approx_vf()
            
            j_local = midas_dict['alpha'] * midas_dict['ug1'] + (1 - midas_dict['alpha']) * midas_dict['vf']
            vgj = midas_dict['ug1'] - j_local
            alpha_j = midas_dict['alpha'] * j_local
            midas_dict.update({'vgj': vgj})
            midas_dict.update({'j': j_local})
            midas_dict.update({'alpha_j': alpha_j})

    return cond.void_area_avg('vgj')

def calc_vr_uncertainty(cond, sigma_vg=0.1, sigma_alpha=0.05, sigma_dp=0.03, percentage = True):
    """Function to calculate the uncertainty in pitot-tube measurements
    
    **Args:**

        - ``sigma_vg``: _description_. Defaults to 0.1.
        - ``sigma_alpha``: _description_. Defaults to 0.05.
        - ``sigma_dp``: _description_. Defaults to 0.03.
        - ``percentage``: are the preceeding . Defaults to True.
    
    Returns:
        - _description_
    """

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            alpha = midas_dict['alpha']
            dp = midas_dict['delta_p']

            if percentage:
                midas_dict['sigma_vg'] = sigma_vg * midas_dict['ug1']
            else:
                midas_dict['sigma_vg'] = sigma_vg

            if percentage:
                midas_dict['sigma_dp'] = sigma_dp * midas_dict['delta_p']
            else:
                midas_dict['sigma_dp'] = sigma_dp

            if percentage:
                midas_dict['sigma_alpha'] = sigma_alpha * midas_dict['alpha']
            else:
                midas_dict['sigma_alpha'] = sigma_alpha
            
            if dp == 0:
                midas_dict['sigma_vf'] = 0
                midas_dict['sigma_vr'] = 0
            else:
                midas_dict['sigma_vf'] = np.sqrt( 1./(2*cond.rho_f) * (sigma_dp**2/((1-alpha)*dp)  + sigma_alpha**2 * dp / (1-alpha)**3) )
                midas_dict['sigma_vr'] = np.sqrt( midas_dict['sigma_vf']**2 + midas_dict['sigma_vg']**2)

    return cond.area_avg('sigma_vr')
