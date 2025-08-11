from .config import *
from .plot_utils import *
from scipy import interpolate
from scipy.optimize import minimize
import warnings

def calc_fric(cond, method = 'Blasius', m = 0.316, n=0.25):
    """Calculates friction factor for each phase based on bulk Re
    
    **Args**:
    
        - ``method``: Method to use to calculate friction factor. Defaults to 'Blasius'. Options:

            - ``'Blasius'``

            .. math:: f_{k} = \\frac{m}{\\text{Re}_{k}^{n}}

        - ``m``: Constant for Blasius-type correlation. Defaults to 0.316.
        - ``n``: Constant for Blasius-type correlation. Defaults to 0.25.
    
    **Raises**:
    
        - ``NotImplementedError``: If invalid method selected
    
    **Returns**:
    
        - Tuple of ``(f_f, f_g)``
        - Stores ``cond.ff``, ``cond.fg``, ``cond.tau_fw``
    """
    
    Re_f = cond.rho_f * cond.jf * cond.Dh / cond.mu_f
    Re_g = cond.rho_g * cond.jgloc * cond.Dh / cond.mu_g

    if method.lower() == 'blasius':
        ff = m / Re_f**n
        fg = m / Re_g**n
    else:
        raise NotImplementedError("Invalid method, try Blasius")

    cond.ff = ff
    cond.fg = fg
    tau_fw = cond.ff/4 * cond.rho_f * cond.jf**2/2
    cond.tau_fw = tau_fw
    return (ff, fg)

def calc_mu_eff(cond, method='Ishii', alpha_peak = 1.0):
    """Method for effective/mixture viscosity
    
    **Args**:
    
        - ``method``: method to use for calculating mixture viscosity. Defaults to 'Ishii'.

            - ``'Ishii'``

            .. math:: \\mu_{eff} = \\mu_{m} = \\mu_{f} (1 - \\frac{\\alpha}{ \\alpha_{max} })^{-2.5 \\alpha_{max} \\frac{\\mu_g + 0.4 \\mu_g}{\\mu_g + \\mu_g} }
            
            - ``'Ishii-AA'``

            .. math:: \\mu_{eff} = \\mu_{m} = \\mu_{f} (1 - \\frac{\\langle \\alpha \\rangle}{ \\alpha_{max} })^{-2.5 \\alpha_{max} \\frac{\\mu_g + 0.4 \\mu_g}{\\mu_g + \\mu_g} }

        - ``alpha_peak``: parameter for some models. Defaults to 1.0.
    
    **Returns**:
    
        - area average effective viscosity
        - stores ``'mu_eff'`` and ``'mu_m'`` in ``midas_dict``
    """

    cond.mirror()
    alpha_avg = cond.area_avg('alpha')
            
    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():

            if method.lower() == 'ishii':
                mu_m = cond.mu_f * (1 - midas_dict['alpha'] / alpha_peak)**(-2.5*alpha_peak * (cond.mu_g + 0.4*cond.mu_f) / (cond.mu_g + cond.mu_f)  )

            elif method.lower() == 'ishii_aa':
                mu_m = cond.mu_f * (1 - alpha_avg / alpha_peak)**(-2.5*alpha_peak * (cond.mu_g + 0.4*cond.mu_f) / (cond.mu_g + cond.mu_f)  )
                
            elif method.lower() == 'avg_void':
                mu_m = cond.mu_f / (1 - alpha_avg)
            else:
                raise(ValueError("Unknown option for calc_mu_eff"))

            if np.real(mu_m) < 0 or np.imag(mu_m) > 0:
                warnings.warn(f"Non-zero or imaginary mu_eff: {angle}, {rstar}, {method}, {midas_dict['alpha']}, {mu_m}. Setting to mu_f")
                mu_eff = cond.mu_f

            mu_eff = mu_m
            mu_m = mu_eff

            midas_dict.update({'mu_eff': mu_eff})
            midas_dict.update({'mu_m': mu_m})
    try:
        return cond.area_avg('mu_eff')
    except:
        return 0

def calc_dpdz(cond, method = 'LM', m = 0.316, n = 0.25, chisholm = 25, k_m = 0.10, L = None, alpha = None, akapower = 0.875):
    """Calculates pressure gradient
    
    **Args**:
    
        - ``method``: Calculation method. Defaults to 'LM'.
            
            - ``'LM'``, Lockhart-martinelli
            -  ``'Kim'``, Kim-modified Lockhart Martinelli
            - ``'akagawa'``, 
        
        - ``m``: option for :any:`calc_fric`. Defaults to 0.316.
        - ``n``: option for :any:`calc_fric`. Defaults to 0.25.
        - ``chisholm``: Chisholm parameter, the C in Lockhart-Martinelli. Defaults to 25.
        - ``k_m``: Minor loss coefficient. Defaults to 0.10.
        - ``L``: Length of restriction, only matters for 'Kim' method. Defaults to None.
        - ``alpha``: void fraction for Akagawa model. Defaults to None.
        - ``akapower``: power for Akagawa model. Defaults to 0.875.
    
    **Raises**:
    
        - ``NotImplementedError``: if unknown option called
    
    **Returns**:
    
        - ``dpdz``
    """
    
    if chisholm == 'Ryan':
        chisholm = 26 - 4.7*np.cos( cond.theta*np.pi/180 )
    
    if method.lower() == 'lm' or method.lower() == 'lockhart' or method.lower() == 'lockhart-martinelli':
        f_f, f_g = cond.calc_fric(m = m, n = n)

        dpdz_f = f_f * 1/cond.Dh * cond.rho_f * cond.jf**2 / 2
        dpdz_g = f_g * 1/cond.Dh * cond.rho_g * cond.jgloc**2 / 2
        chi2 = dpdz_f / dpdz_g

        phi_f2 = 1 + chisholm/np.sqrt(chi2) + 1 / chi2
        dpdz = phi_f2 * dpdz_f

    elif method.lower() == 'kim':
        f_f, f_g = cond.calc_fric(m = m, n = n)

        dpdz_f = f_f * 1/cond.Dh * cond.rho_f * cond.jf**2 / 2
        dpdz_g = f_g * 1/cond.Dh * cond.rho_g * cond.jgloc**2 / 2
        dpdz_m = k_m * cond.rho_f * cond.jf**2 / 2 / L

        chi2 = dpdz_f / dpdz_g
        chiM2 = dpdz_f / dpdz_m

        phi_f2 = (1 + 1 / chiM2) + np.sqrt(1 + 1 / chiM2) * chisholm / np.sqrt(chi2) + 1 / chi2
        dpdz = phi_f2 * dpdz_f

    elif method.lower() == 'akagawa':
        f_f, f_g = cond.calc_fric(m = m, n = n)

        dpdz_f = f_f * 1/cond.Dh * cond.rho_f * cond.jf**2 / 2
        phi_f2 = ((1 - alpha)**(-akapower))**2

        dpdz = phi_f2 * dpdz_f

    else:
        raise NotImplementedError(f'{method} is not a valid option for calc_dpdz. Try "LM" ')
    
    cond.dpdz = dpdz
    cond.tau_w = dpdz * cond.Dh / 4

    return dpdz

def calc_cd(cond, method='Ishii-Zuber', vr_cheat = False, limit = 10**-6, const_CD = 0.44):
    """Method for calculating drag coefficient
    
    **Args**:
    
        - ``method``: what method to use for modeling :math:`C_{D}`. Defaults to 'Ishii-Zuber'. Options:
            
            - ``'Ishii-Zuber'``

            .. math:: C_{D} = \\frac{24}{\\text{Re}_{b}} (1 + 0.1 \\text{Re}_{b}^{0.75})

            - ``'Schiller-Naumann'``

            .. math:: C_{D} = \\frac{24}{\\text{Re}_{b}} (1 + 0.15 \\text{Re}_{b}^{0.687})

            - ``'constant'`` or ``'const'``

        - ``vr_cheat``: flag to use "vr" from ``midas_dict`` or "vr_model" when calculating :math:`Re_{b}`. Defaults to False.
        - ``limit``: Option to specify a minimum :math:`C_{D}`. Will still calculate :math:`C_{D}` based on ``method``, but the result won't be below this minimum. Defaults to 10**-6. Options:

            - ``'tomiyama'``, 

            .. math:: C_{D, min} = \\frac{8}{3} \\frac{\\text{Eo}}{\\text{Eo} + 4}

            - ``'ishii-chawla'``

            .. math:: C_{D, min} = \\min{(\\frac{2}{3} \\sqrt{\\text{Eo}}, \\frac{8}{3})}

            - User-specified constant

        - ``const_CD``: Constant drag coefficient, if ``method = 'constant'`` . Defaults to 0.44.
    
    **Raises**:
    
        - ``NotImplementedError``: If unknown method called
    
    **Returns**:
    
        - area-averaged drag coefficient
        - Stores ``'cd'`` and ``'Reb'`` in ``midas_dict``
    """

    cond.calc_mu_eff()

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():

            if vr_cheat:
                Reb = (1 - midas_dict['alpha']) * midas_dict['Dsm1'] * cond.rho_f * abs(midas_dict['vr']) / midas_dict['mu_m']
            else:

                if 'vr_model' not in midas_dict.keys(): # Initialize for iteration
                    midas_dict.update({'vr_model': -1})

                Reb = (1 - midas_dict['alpha']) * midas_dict['Dsm1'] * cond.rho_f * abs(midas_dict['vr_model']) / midas_dict['mu_m']

            try:
                if Reb < 0:
                    warnings.warn("RuntimeWarning: Reb negative!")
                    print(Reb, midas_dict['alpha'], midas_dict['Dsm1'], abs(midas_dict['vr_model']), midas_dict['mu_m'])
                else:
                    pass
            except:
                warnings.warn("RuntimeWarning: Reb imaginary!")
                print(Reb, midas_dict['alpha'], midas_dict['Dsm1'], abs(midas_dict['vr_model']), midas_dict['mu_m'])
            
            midas_dict.update({'Reb': Reb})

            if method == 'Ishii-Zuber' or method == 'IZ' or method == 'Ishii':

                if Reb > 0:
                    cd = 24/Reb * (1 + 0.1*Reb**0.75)
                else:
                    cd = limit

            elif method.lower() == 'schiller-naumann':
                Reb = (1 - midas_dict['alpha']) * midas_dict['Dsm1'] * cond.rho_f * abs(midas_dict['vr_model']) / cond.mu_f

                cd = 24/Reb * (1 + 0.15*Reb**0.687)
            
            elif method == 'constant' or method == 'const':
                cd = const_CD

            if type(limit) == str:
                if limit.lower() == "tomiyama":
                    eo = cond.g * (cond.rho_f - cond.rho_g) * midas_dict['Dsm2']
                    limit = 8/3 * eo / (eo + 4)
                    midas_dict.update({'eo': eo})
                    
                elif limit.lower() == 'ishii-chawla':
                    eo = cond.g * (cond.rho_f - cond.rho_g) * midas_dict['Dsm2']
                    limit = min(2/3*np.sqrt(eo), 8/3)
                    midas_dict.update({'eo': eo})
                
                else:
                    raise NotImplementedError(f"{limit} not a valid type for limiting behavior. Please enter tomiyama, Ishii-Chawla, or set a constant limit (e.g. limit = 0.44)")
            
            cd = max(limit, cd) # Either 0, set by user, or set by above string

            midas_dict.update({'cd': cd})

    return cond.area_avg('cd')

def calc_cl(cond, method='tomiyama', sharma_factor = False):
    """TODO, not implemented
    
    **Args**:
    
        - ``method``: method to use. Defaults to 'tomiyama'.

            - ``'tomiyama'``
            - ``'hibiki-ishii'``

        - ``sharma_factor``: not implemented. Defaults to False.
    """
    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            
            if method.lower() == 'tomiyama':
                CL = 0
            elif method.lower() == 'hibiki-ishii' or method.lower() == 'hi':
                CL = 0

            midas_dict.update({'CL': CL})

def calc_vr_model(cond, method='km1_simp', kw = 0.654, n=1, Lw = 5, kf = 0.113, 
                    iterate_cd = True, initial_vr = None, 
                    quiet = True, recalc_cd = True, iter_tol = 1e-4, custom_f = None, CC = 1):
    """_summary_
    
    **Args**:
    
        - ``method``: what method to use for modeling :math:`v_{r}`. Defaults to 'km1_simp'. Options include:

            - ``'wake_1'``, depracated, earliest attempt, don't use

            .. math:: v_{r} = K_{w} v_{f} C_{D}^{1/3}

            - ``'wake_alpha'``, depracated, don't use

            .. math:: v_{r} = K_{w} (1-\\alpha)^{n} v_{f} C_{D}^{1/3}

            - ``'wake_alpha2'``, depracated, don't use

            .. math:: v_{r} = K_{w} \\alpha (1-\\alpha)^{n} v_{f} C_{D}^{1/3}

            - ``'wake_lambda'``, depracated, don't use
            
            .. math:: v_{r} = K_{w} v_{f} C_{D}^{1/3} \\left(\\frac{D_{sm1}}{\\lambda}\\right)^{2/3}

            - ``'wake_vg_lambda'``, depracated, don't use
            
            .. math:: v_{r} = K_{w} \\alpha (1-\\alpha)^{n} v_{f} C_{D}^{1/3}

            - ``'wake_vr'`` or ``'hubris'``, depracated, don't use. Ugly integral
            - ``'km1'``, depracated, don't use
            
            .. math:: v_{r} = K_{w} \\frac{\\pi}{4} \\alpha v_{f} C_{D}^{1/3} \\frac{2^{-1/3}-L_{w}^{1/3}}{0.5-L_{w}} + K_{f} v_{f}

            - ``'prelim'`` or ``'km1_simp'`` or ``'final_horizontal'``, model used for adix prelim. Recommended for horizontal flows
            
            .. math:: v_{r} = - K_{w} \\alpha v_{f} C_{D}^{1/3} - K_{f} v_{f}

            - ``'prelim_plus'`` , model with gravity and friction
            
            .. math:: v_{r} = K_{w} \\alpha v_{f} C_{D}^{1/3} - K_{f} v_{f} + \\sqrt{ \\frac{4}{3} \\frac{ \\langle \\langle D_{sm1} \\rangle \\rangle }{C_{D}} \\left( f_{f} \\frac{j_{f}^{2}}{2} + (1-\\alpha) \\left( 1-\\frac{\\rho_{g}}{\\rho_{f}} \\right) g_{z} \\right)}

            - ``'final'``, final model from adix, with gravity. Recommended for horizontal flows
            
            .. math:: v_{r} = \\text{sgn} \\left(v_{r}\\right ) K_{w} \\alpha v_{f} C_{D}^{1/3} - K_{f} v_{f} + \\sqrt{ \\frac{4}{3} \\frac{ \\langle \\langle D_{sm1} \\rangle \\rangle }{C_{D}} \\left( (1-\\alpha) \\left( 1-\\frac{\\rho_{g}}{\\rho_{f}} \\right) g_{z} \\right)}

            - ``'ishii-chawla'``

            .. math:: v_{r} = \\sqrt{2} \\left ( \\frac{\\sigma g (\\rho_{f} - \\rho_{g})}{\\rho_{f}^{2}} \\right)^{0.25}
            
            - ``'proper_integral'``, depracated, do not use
            - ``'proper_integral_alpha'``, depracated, do not use
            - ``'Chahed'``, see Chahed et al. (2017), not recommended

            .. math:: v_{r} |v_{r}| = \\frac{4}{3} \\frac{D_{sm,1}}{C_{D}} \\left ( \\frac{4 \\tau_fw}{\\rho_{f}} + g_{z} (1-\\alpha) - \\frac{C_{C}}{\\alpha} \\frac{\\partial \\alpha}{\\partial r} \\frac{\\partial v_{f}}{\\partial r} \\right)

        - ``kw``: wake coefficient. Defaults to 0.654.
        - ``n``: power. Defaults to 1.
        - ``Lw``: Effective wake length. Defaults to 5.
        - ``kf``: Liquid coefficient. Defaults to 0.113.
        - ``iterate_cd``: Option to iterate :math:`C_{D}` with the newly calculated :math:`v_{r}`. Defaults to True.
        - ``initial_vr``: Initialize :math:`v_{r}` for iteration. Defaults to None.
        - ``quiet``: output messages. Defaults to True.
        - ``recalc_cd``: Recalculate :math:`C_{D}`. Defaults to True.
        - ``iter_tol``: Convergence tolerance for iteration. Defaults to 1e-4.
        - ``custom_f``: Custom friction factor, if using a method that uses the friction factor. Defaults to None.
        - ``CC``: Chahed constant. Defaults to 1.
    
    **Raises**:
    
        - ``ValueError``: if invalid method selected.
    
    **Returns**:
    
        - area average relative velocity calculated by the model
    """
    

    MAX_ITERATIONS = 100
    iterations = 0
    initialize_vr = True

    if initial_vr is None:
        if cond.theta == 90:
            initial_vr = 0.5
        else:
            initial_vr = -0.5

    while True:
        if recalc_cd:
            if iterate_cd:
                
                if initialize_vr:
                    for angle, r_dict in cond.data.items():
                        for rstar, midas_dict in r_dict.items():
                            midas_dict.update(
                                {'vr_model': initial_vr}
                            )
                    initialize_vr = False

                cond.calc_cd(vr_cheat=False)
            else:
                cond.calc_cd(vr_cheat=True)
        try:
            old_vr = cond.area_avg('vr_model', recalc=True)
        except KeyError:
            old_vr = initial_vr # Initialize?

        vr_name = "vr_" + method

        for angle, r_dict in cond.data.items():
            for rstar, midas_dict in r_dict.items():
                
                if method == 'wake_1':
                    vr = kw  * midas_dict['vf'] * midas_dict['cd']**(1./3)
                
                elif method == 'wake_alpha':
                    vr = kw  * (1 - midas_dict['alpha'])**n * midas_dict['vf'] * midas_dict['cd']**(1./3)
                
                elif method == 'wake_alpha2':
                    vr = kw  * (midas_dict['alpha']*(1 - midas_dict['alpha']))**n * midas_dict['vf'] * midas_dict['cd']**(1./3)

                elif method == 'wake_lambda':
                    cond.calc_avg_lat_sep()
                    vr = kw  * midas_dict['vf'] * midas_dict['cd']**(1./3) * midas_dict['Dsm1']**(2/3) * midas_dict['lambda']**(-2./3)

                elif method == 'wake_vg_lambda':
                    cond.calc_avg_lat_sep()
                    vr = midas_dict['ug1'] /(1 + kw * midas_dict['cd']**(1./3) * midas_dict['Dsm1']**(2/3) * midas_dict['lambda']**(-2./3) )

                elif method == 'hubris' or method == 'wake_vr':
                    cd = midas_dict['cd']
                    c2 = kw
                    vr = 1/2* midas_dict['ug1']*(1 - 3 *kw* cd**(1/3)*(2**(2/3) - 2* Lw**(1/3)) - 2* Lw + 6* c2**(3/2)* np.sqrt(cd) *(np.arctan(1./(2**(1/3) * np.sqrt(c2) * cd**(1/6))) - np.arctan(1/((np.sqrt(c2)* cd**(1/6))/Lw**(1/3)))))
                
                elif method == 'km1':
                    vr = kw * (np.pi/4)**(1/3) * midas_dict['alpha'] * midas_dict['vf'] * midas_dict['cd']**(1./3) *  (2**(-1./3) - Lw**(1/3))/(0.5 - Lw) + kf * midas_dict['vf']
                
                elif method == 'km1_simp' or method == 'prelim' or method.lower() == 'final_horizontal':
                    vr = -kw * midas_dict['alpha'] * midas_dict['vf'] * midas_dict['cd']**(1./3) - kf * midas_dict['vf']

                elif method == 'prelim_plus':
                    if custom_f is None:
                        ff, fg = cond.calc_fric()
                    else:
                        ff = custom_f

                    # if midas_dict['Dsm1'] == 0 and rstar != 1.0:
                    #     print(cond, angle, rstar)

                    try:
                        vr = (
                        +kw * midas_dict['alpha'] * midas_dict['vf'] * midas_dict['cd']**(1./3) - kf * midas_dict['vf'] 
                        + np.sqrt( 4./3 * cond.void_area_avg('Dsm1')*0.001/midas_dict['cd'] * ( ff/cond.Dh * cond.jf**2/2 + 
                                                                                (1 - midas_dict['alpha'])*(1-cond.rho_g/cond.rho_f) * cond.gz ) )
                        )
                    except ZeroDivisionError:
                        vr = 0

                    if vr == np.inf and midas_dict['cd'] == 0:
                        if not quiet: warnings.warn(f"vr nan for {angle, rstar}, setting to 0")
                        vr = 0
                    
                    if vr != vr:
                        if not quiet: warnings.warn(f"vr nan for {angle, rstar}, setting to 0")
                        vr = 0

                elif method == 'final':
                    try:
                        vr = (
                        np.sign(old_vr) * kw * midas_dict['alpha'] * midas_dict['vf'] * midas_dict['cd']**(1./3) - kf * midas_dict['vf'] 
                        + np.sqrt( 4./3 * cond.void_area_avg('Dsm1')*0.001/midas_dict['cd'] * (1 - midas_dict['alpha'])*(1-cond.rho_g/cond.rho_f) * cond.gz ) )
                    except ZeroDivisionError:
                        vr = 0
                    

                elif method == 'proper_integral':
                    warnings.warn("This method is probably no good, messed up the math")
                    vr = midas_dict['ug1'] / ( (0.5 - Lw) + kw * midas_dict['cd']**(1./3) *(np.pi/4)**(1/3)* (2**(-1./3) - Lw**(1/3)))

                elif method == 'proper_integral_alpha':
                    warnings.warn("This method is probably no good, messed up the math")
                    vr = midas_dict['ug1'] / ( (0.5 - Lw) + kw * midas_dict['alpha']**n *midas_dict['cd']**(1./3) *(np.pi/4)**(1/3)* (0.5**(1./3) - Lw**(1/3)))

                elif method.lower() == 'ishii-chawla' or method.lower() == 'ishii' or method.lower() == 'ishii chawla':
                    vr = np.sqrt(2) * (cond.sigma * cond.g * (cond.rho_f - cond.rho_g) / (cond.rho_f**2))**0.25
                
                elif method.lower() == 'chahed': # see Chahed et al. (2017)
                    cond.calc_dpdz()
                    cond.calc_grad('alpha')
                    cond.calc_grad('vf')
                    if abs(midas_dict['alpha']) < 1e-6 or abs(midas_dict['cd']) < 1e-6:
                        discrim = 0
                    else:
                        # print(cond.Dh, cond.rho_f, midas_dict['cd'], midas_dict['alpha'])
                        discrim = 4/3 * midas_dict['Dsm1'] / midas_dict['cd'] * ( 4/cond.Dh * cond.tau_w  / cond.rho_f + cond.gz*(1-midas_dict['alpha']) - CC/midas_dict['alpha'] * midas_dict['grad_alpha_r'] * midas_dict['grad_vf_r'])

                    vr = float(np.sign(discrim) * np.sqrt(discrim))
                else:
                    print(f"{method} not implemented")
                    return -1

                if rstar == 1:
                    vr = 0

                if vr > 2*cond.jf:
                    vr = 2*cond.jf
                elif vr < -2*cond.jf:
                    vr = -2*cond.jf
        
                midas_dict[vr_name] = vr
                midas_dict['vr_model'] = vr


        iterations += 1

        try:
            cond.area_avg('vr_model', recalc=True)
        except:
            print(f"Error calculating vr" )

        if old_vr == np.inf or old_vr != old_vr or vr == np.inf:
            print(f"Error calculating vr: {cond.area_avg('vr_model', recalc=True)}")

        if old_vr == 0:
            if not quiet:
                print(f"vr_model calculated as 0 after {iterations} iterations")
                print(old_vr, cond.area_avg('vr_model', recalc=True))
            return

        if abs(old_vr - cond.area_avg('vr_model', recalc=True)) / abs(old_vr) < iter_tol:
            if not quiet:
                print(f"vr_model converged in {iterations} iterations")
                print(old_vr, cond.area_avg('vr_model', recalc=True))
            return cond.area_avg("vr_model")
        
        if iterations > MAX_ITERATIONS:
            print("Warning, max iterations exceeded in calculating vr_model")
            print(f"{old_vr}\t{cond.area_avg('vr_model', recalc=True)}\t{(old_vr - cond.area_avg('vr_model', recalc=True))/old_vr*100}")
            return
        
    
    return cond.area_avg("vr_model")

def calc_vgj_model(cond):
    """Method for calculating local Vgj based on models

    Stores:
        - ``midas_dict['vgj_model'] = (1 - midas_dict['alpha']) * midas_dict['vr_model']``
    
    """

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():

            if 'vr_model' not in midas_dict.keys():
                print(f"Warning: vr_model not found for {angle}, {rstar}, calling calc_vr_model with default inputs")
                cond.calc_vr_model()

            midas_dict['vgj_model'] = (1 - midas_dict['alpha']) * midas_dict['vr_model']

    return

def calc_IS_term(cond, method = 'power', n=2, mu = 1.5):
    """Calculate the interfacial shear term, for the 1-D averaged TFM, :math:`\\nabla \\alpha \\bullet \\tau_{i}`
    
    **Args**:

        - ``method``: method used to calculate term. Defaults to 'power'. Options:

            - ``power``

            .. math:: \\tau_{i} = \\tau_{fw} (r^{*})^{n}

            - ``power_total``
            - ``lognorm``
            - ``alpha``
        - ``n``: power. Defaults to 2.
        - ``mu``: mean of lognormal distribution. Defaults to 1.5.
    
    **Returns**:

        - ``cond.area_avg('ISxgrad')``. ``'ISxgrad`` and ``'IS'`` stored in ``midas_dict``
    """
    cond.calc_cd()
    cond.calc_fric()
    cond.calc_grad('alpha')
    
    for angle, r_dict in cond.data.items():
        for r_star, midas_dict in r_dict.items():
                
                if r_star == 0 and n < 0:
                    midas_dict.update(
                        {'IS': 0 , 'ISxgrad': 0}
                        )
                    continue

                if method == 'power':
                    taui = float(cond.tau_fw * r_star**n)
                    midas_dict.update(
                        {'ISxgrad': float(taui * ( midas_dict['grad_alpha_r'] )) , 'IS': taui}
                        )
    
                elif method == 'power_total':
                    taui = float(cond.tau_fw * r_star**n)
                    midas_dict.update(
                        {'ISxgrad': float(taui * (midas_dict['grad_alpha_total'] )) , 'IS': taui}
                    )
                    
                elif method == 'lognorm':
                    taui = cond.tau_fw * np.exp( -(np.log(1-r_star) + mu)**2 )
                    midas_dict.update(
                        {'ISxgrad': float(taui * (midas_dict['grad_alpha_r'] )) , 'IS': float(taui)}
                    )

                elif method == 'alpha':
                    taui = n * cond.tau_fw * midas_dict['alpha'] / cond.area_avg('alpha')
                    midas_dict.update(
                        {'ISxgrad': float(taui * (midas_dict['grad_alpha_r'] )) , 'IS': float(taui)}
                    )


    return cond.area_avg('ISxgrad')

def calc_aa_vr_model(cond, method='km1_naive', IS_method = 'power', kw=0.654, kf=0.113, Lw = 5, Ctau=1, n=2, IS_mu = 1.5, Cvfacd = 1):
    """to calculate estimate the area-averaged relative velocity
    
    **Args:**

        - ``method``: method to use. Defaults to 'km1_naive'.

            - ``km1_naive``, depracated, similar to ``prelim`` but with an older form and different coefficients

            .. math:: \\langle v_{r} \\rangle = K_{f} \\frac{\\langle j_{f} \\rangle}{1 - \\langle \\alpha \\rangle} + K_{w} \\frac{\\pi}{4} \\frac{2^{-1/3}-L_{w}^{1/3}}{0.5-L_{w}} \\langle C_{D} \\rangle^{1/3} \\langle \\alpha \\rangle \\frac{\\langle j_{f} \\rangle}{1 - \\langle \\alpha \\rangle}

            - ``km1_naive2``, depracated, same as ``prelim``, but signs of :math:`K_{w}` and :math:`K_{f}` flipped.
            - ``prelim``, Same as final, but no covariance. For backwards compatability 

            .. math:: \\langle v_{r} \\rangle = - K_{f} \\frac{\\langle j_{f} \\rangle}{1 - \\langle \\alpha \\rangle} - K_{w} \\langle C_{D} \\rangle^{1/3} \\langle \\alpha \\rangle \\frac{\\langle j_{f} \\rangle}{1 - \\langle \\alpha \\rangle}

            - ``final``, recommended method, Eq. (4.37) of Dix (2025)

            .. math:: \\langle v_{r} \\rangle = - K_{f} \\frac{\\langle j_{f} \\rangle}{1 - \\langle \\alpha \\rangle} - K_{w} C_{\\alpha, v_{f}} \\langle C_{D} \\rangle^{1/3} \\langle \\alpha \\rangle \\frac{\\langle j_{f} \\rangle}{1 - \\langle \\alpha \\rangle}

            - ``IS_Ctau``, uses a 1-D interfacial shear model with :math:`C_{\\tau}` to calculate :math:`v_{r}`. Not recommended

            .. math:: \\langle v_{r} \\rangle | \\langle v_{r} \\rangle| = \\frac{8 r_{b}}{3} \\frac{1}{C_{D} \\rho_{f}} ((1-C_{\\tau})\\frac{4 \\tau_{fw}}{D_{h}} + (1-\\langle \\alpha \\rangle)g_{z}(\\rho_{f} - \\rho_{g}) )
            
            - ``IS``, using 1-D interfacial shear model. Uses :meth:`~MARIGOLD.Condition.Condition.calc_IS_term` for :math:`\\tau_{i}`. Not recommended

            .. math:: \\langle v_{r} \\rangle | \\langle v_{r} \\rangle| = \\frac{8 r_{b}}{3} \\frac{1}{C_{D} \\rho_{f}} (\\frac{4 \\tau_{fw}}{D_{h}} + (1-\\langle \\alpha \\rangle)g_{z}(\\rho_{f} - \\rho_{g}) - \\frac{1}{\\langle \\alpha \\rangle} \\tau_{i})
        
        - ``IS_method``: When ``method='IS'``, what method to pass to use for interfacial shear calculation. See :meth:`~MARIGOLD.Condition.Condition.calc_IS_term` for options. Defaults to ``'power'``. 
        - ``kw``: Wake coefficient. Defaults to 0.654.
        - ``kf``: Liquid coefficient. Defaults to 0.113.
        - ``Lw``: Effective wake length. Defaults to 5.
        - ``Ctau``: Interfacial shear constant. Defaults to 1.
        - ``n``: Power constant. Defaults to 2.
        - ``IS_mu``: Interfacial shear viscosity. Defaults to 1.5.
        - ``Cvfacd``: Covariance. Defaults to 1.
    
    Returns:
        - ``cond.aa_vr``
    """

    if method == 'km1_naive':
        vr = kw * (np.pi/4)**(1/3) * cond.area_avg('alpha') * cond.jf / (1 - cond.area_avg('alpha')) * cond.area_avg('cd')**(1./3) *  (2**(-1./3) - Lw**(1/3))/(0.5 - Lw) + kf * cond.jf / (1 - cond.area_avg('alpha'))
    
    elif method == 'km1_naive2' :#or method == 'prelim':
        cond.calc_cd()
        vr = kw * cond.area_avg('alpha') * cond.jf / (1 - cond.area_avg('alpha')) * cond.area_avg('cd')**(1./3)  + kf * cond.jf / (1 - cond.area_avg('alpha'))

    elif method == 'km1_simp' or method == 'prelim':
        vr = -kw * cond.area_avg('alpha') * cond.jf / (1 - cond.area_avg('alpha')) * cond.area_avg('cd')**(1./3) - kf * cond.jf / (1 - cond.area_avg('alpha'))
    
    elif method == 'final':
        cond.calc_cd()
        vr = Cvfacd * -kw * cond.area_avg('alpha') * cond.jf / (1 - cond.area_avg('alpha')) * cond.area_avg('cd')**(1./3) - kf * cond.jf / (1 - cond.area_avg('alpha'))
    
    elif method == 'IS_Ctau':
        cond.calc_cd()
        cond.calc_fric()
        rb = cond.void_area_avg('Dsm1') / 2 /1000 # Convert to m
        CD = cond.void_area_avg('cd')
        alpha = cond.area_avg('alpha')

        discrim = (1-alpha)*cond.gz * (cond.rho_f - cond.rho_g) + (1-Ctau)*4*cond.tau_fw/cond.Dh
        vr = np.sign(discrim) * np.sqrt(8*rb/3 * 1/(CD * cond.rho_f) * abs( discrim ))

    elif method == 'IS':
        IS_term = cond.calc_IS_term(method = IS_method, n = n, mu = IS_mu)

        rb = cond.void_area_avg('Dsm1') / 2 /1000 # Convert to m
        CD = cond.void_area_avg('cd')
        alpha = cond.area_avg('alpha')

        discrim = 4*cond.tau_fw/cond.Dh + (1-alpha)*cond.gz * (cond.rho_f - cond.rho_g) - 1/alpha * IS_term
        vr = np.sign(discrim) * np.sqrt(8*rb/3 * 1/(CD * cond.rho_f) * abs( discrim ))
    else:
        raise(ValueError("Invalid calc_aa_vr_model method"))
        

    cond.vwvgj = (1-cond.area_avg('alpha'))*vr # Legacy, do not use
    cond.aa_vr = vr
    return vr

def calc_vw_aa_Vgj_model(cond, Kw=0.654, Kf=0.113):
    """ Calculate drift velocity, :math:`\\langle \\langle V_{gj} \\rangle \\rangle`, based on Dix (2025) model. See Eq. (4.47) of his thesis
    
    **Args:**

        - ``Kw``: Wake coefficient. Defaults to 0.654.
        - ``Kf``: Liquid coefficient. Defaults to 0.113.
    
    Returns:
        - ``cond.vw_aa_Vgj_model``
    """
    cond.Cajf = 0.12 * cond.jf + 0.43* cond.jgloc**-0.1
    cond.Cajc = 61.4 * cond.jf**-2.25

    for angle, rdict in cond.data.items():
        for rstar, midas_dict in rdict.items():
            midas_dict.update({'cd13': midas_dict['cd']**(1./3)})

    cond.vw_aa_Vgj_model = -Kf * cond.Cajf * cond.jf - Kw * (cond.Cajc * cond.Cajf * cond.area_avg('alpha') * cond.jf * cond.void_area_avg('cd13'))
    
    return cond.vw_aa_Vgj_model
