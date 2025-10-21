from .config import *

def avg(cond, param: str, include_zero=True) -> float:
    """Calculates a basic average of a parameter

    This function does not mirror the data. Note that mirroring might change the average.
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter to average. See :func:`~MARIGOLD.Condition.print_params` for options
        - ``include_zero``: include zero values . Defaults to True.
    
    **Returns**:
    
        - the average
    """

    count = 0
    avg_param = 0
    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            
            if include_zero:
                count += 1
                avg_param += midas_dict[param]
            else:
                if abs(midas_dict[param]) > 0:
                    count += 1
                    avg_param += midas_dict[param]

    return avg_param / count

def area_avg(cond, param: str, even_opt='first', recalc=True, method=None) -> float:
    """Method for calculating the area-average of a parameter

    .. math:: \\langle \\psi \\rangle = \\frac{1}{A} \iint_A \\psi(r,\\varphi) r \,dr\,d\\varphi
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter to area-average. See :func:`~MARIGOLD.Condition.print_params` for options
        - ``even_opt``: option for ``scipy.integrate.simpsons``. Defaults to 'first'.
        - ``recalc``: if true, will recalculate area average. Defaults to True.
        - ``method``: method to area-average. Defaults to None.

            - ``legacy``, using the same method as the Excel spreadsheets
            - ``legacy_old``, actually what we use in spreadsheets, hardcoded values
            - None, will use ``scipy.integrate.simpsons``. Recommended option
    
    **Raises**:
    
        - ``KeyError``: if ``param`` not found
    
    **Returns**:
    
        - area-averaged value
    """

    if not cond.check_param(param):
        raise KeyError(f"Invalid parameter {param} selected. Not present at {cond.check_param_loc(param)}")
    
    if (param in cond.area_avgs.keys()) and (not recalc):
        return cond.area_avgs[param] # why waste time, if we already calculated this don't do it again
    
    if method == 'legacy':
        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        if not cond.mirrored:
            warnings.warn("Mirroring in area-avg")
            cond.mirror()

        for angle, r_dict in cond.data.items():

            if angle == 360:    # We already have 0
                continue

            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi / 180)
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    try:
                        rs_temp.append(rstar)
                        vars_temp.append(float(midas_dict[param] * rstar))
                    except:
                        if debug: print('Problem with:', angle, r_dict, param)
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            vars.reverse()      # I noticed too late that the list goes from r/R = 0 to 1, not the other way around
            rs.reverse()        # Already wrote the actual Simpson's rule part, easier to just do this

            if debug: print("Arrays to integrate", angle, rs, vars, file=debugFID)

            if len(rs) != len(vars):
                ValueError( f"rs to integrate over {rs} must be the same length as params {vars}, occured at {angle}" )

            delta = abs(np.diff(rs))                # r/R steps

            la = 0
            for idx, var in enumerate(vars):
                coeff = 2 * (2**(idx % 2)) / 3      # Simpson's Rule coefficient
                
                if idx < len(delta):
                    if idx > 0:
                        if delta[idx - 1] == delta[idx]:
                            S = delta[idx] * coeff * var
                        else:
                            coeff = coeff / 2
                            S = (delta[idx - 1] * coeff * var) + (delta[idx] * coeff * var)
                    else:
                        S = delta[idx] * coeff * var
                else:
                    S = delta[idx - 1] * coeff * var
                
                la = la + S
            
            param_r.append(la)

        I = sum(param_r) / 8
        cond.area_avgs.update({param: I})

    elif method == 'legacy_old':
        # There's definitely a more elegant way to do this, but I just want Talley's data to work for now

        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        if not cond.mirrored:
            warnings.warn("Mirroring in area-avg")
            cond.mirror()

        for angle, r_dict in cond.data.items():
            if angle == 360:
                continue

            if 0.95 in r_dict:
                S_1 = 0.05 * sum((1 * 0,
                        4 * 0.95 * r_dict[0.95][param],
                        2 * 0.90 * r_dict[0.90][param],
                        4 * 0.85 * r_dict[0.85][param],
                        1 * 0.80 * r_dict[0.80][param],
                        )) / 3
                
                S_2 = 0.10 * sum((1 * 0.80 * r_dict[0.8][param],
                        4 * 0.70 * r_dict[0.7][param],
                        2 * 0.60 * r_dict[0.6][param],
                        4 * 0.50 * r_dict[0.5][param],
                        1 * 0.40 * r_dict[0.4][param],
                        )) / 3
                
                S_3 = 0.20 * sum((1 * 0.40 * r_dict[0.4][param],
                        4 * 0.20 * r_dict[0.2][param],
                        2 * 0.00 * r_dict[0.0][param],
                        )) / 3
            else:
                S_1 = 0

                S_2 = 0.10 * sum((1 * 0,
                        4 * 0.90 * r_dict[0.9][param],
                        2 * 0.80 * r_dict[0.8][param],
                        4 * 0.70 * r_dict[0.7][param],
                        2 * 0.60 * r_dict[0.6][param],
                        4 * 0.50 * r_dict[0.5][param],
                        1 * 0.40 * r_dict[0.4][param],
                        )) / 3
                
                S_3 = 0.20 * sum((1 * 0.40 * r_dict[0.4][param],
                        4 * 0.20 * r_dict[0.2][param],
                        2 * 0.00 * r_dict[0.0][param],             # Might be doubling up
                        )) / 3

            param_r.append(sum((S_1, S_2, S_3)))

        I = sum(param_r) / 8
        cond.area_avgs.update({param: I})

    else:
        # We have to integrate twice, once with resepect to r, again with respect to phi
        # Start with r
        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        if not cond.mirrored:
            warnings.warn("Mirroring in area-avg")
            cond.mirror()

        for angle, r_dict in cond.data.items():

            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0: # This should be unnecessary now with the new mirror, but it's not hurting anyone by being here
                    try:
                        rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                        vars_temp.append( float( midas_dict[param] * rstar)) # Floatify to avoid np inhomogeneous array issues
                    except:
                        if debug: print('Problem with:', angle, r_dict, param)
                    #if debug: print(angle, midas_dict, file=debugFID)
            
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            if debug: print("Arrays to integrate", angle, rs, vars, file=debugFID)

            if len(rs) != len(vars):
                ValueError( f"rs to integrate over {rs} must be the same length as params {vars}, occured at {angle}" )
                
            try:
                param_r.append( integrate.simpson(y=vars, x=rs) ) # Integrate wrt r
            except Exception as e:
                if debug:
                    print(e)
                    print(rs, vars)
            if debug: print("calculated integral:", integrate.simpson(y=vars, x=rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)

        param_r = [param for _, param in sorted(zip(angles, param_r))]
        angles = sorted(angles)

        I = integrate.simpson(y=param_r, x=angles) / np.pi # Integrate wrt theta, divide by normalized area
        cond.area_avgs.update({param: I})

    return I

def circ_segment_area_avg(cond, param:str, hstar:float, int_err = 10**-4) -> float:
    """Area-average over a circular segment defined by :math:`h^{*}`

    *WARNING*: uses ``scipy.integrate.dblquad``. May be computationally expensive
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        - ``hstar``: :math:`h^{*}` defining the circular segment of interest. See :any:`find_hstar_pos`
        - ``int_err``: acceptable warning for integral error (passed to scipy.integrate.dblquad). Defaults to 10**-4.
    
    **Raises**:
    
        - ``KeyError``: If ``param`` not available
    
    **Returns**:
    
        - area-average over the circular segment
    """

    # Check that the parameter that the user requested exists
    if not cond.check_param(param):
        raise KeyError(f"Invalid parameter {param} selected. Not present at {cond.check_param_loc(param)}")
    
    cond.mirror()

    rs = []
    phis = []
    vals = []
    for angle, r_dict in cond.data.items():
        for r_star, midas_dict in r_dict.items():
            rs.append(r_star)
            phis.append(angle * np.pi/180) # Convert degrees to radians
            try:
                vals.append(midas_dict[param] * r_star) # Don't forget r dr dθ
            except:
                vals.append(np.NaN * r_star) # Don't forget r dr dθ
                print(f"Could not find {param} for φ = {angle}, r = {r_star}. Substituting NaN")

    # Set up interpolation
    rs = np.asarray(rs)
    phis = np.asarray(phis)
    vals = np.asarray(vals)   

    interp = interpolate.LinearNDInterpolator(list(zip(rs, phis)), vals) # TODO refactor to use __call__

    # To integrate, need to get the bounds. They're different depending on if h* > 1
    def lower_r_bound_pos(phi):
        return max((1 - hstar) / np.sin(phi), 0)

    def upper_r_bound_neg(phi):
        if (phi <= 3*np.pi/2 - np.arccos(hstar-1)) or (phi >= 3*np.pi/2 + np.arccos(hstar-1)):
            return 1
        else:
            return (1 - hstar) / np.sin(phi)

    def dumb_func(r, phi): # For area integration
        return r

    if hstar <= 1:
        I = integrate.dblquad(interp, np.pi/2 - np.arccos(1-hstar), np.pi/2 + np.arccos(1-hstar), lower_r_bound_pos, 1, epsabs=int_err )[0] / integrate.dblquad(dumb_func, np.pi/2 - np.arccos(1-hstar), np.pi/2 + np.arccos(1-hstar), lower_r_bound_pos, 1, epsabs=int_err )[0]
    elif hstar > 1:
        I = integrate.dblquad(interp, 0, 2*np.pi, 0, upper_r_bound_neg, epsabs=int_err )[0] / integrate.dblquad(dumb_func, 0, 2*np.pi, 0, upper_r_bound_neg, epsabs=int_err )[0]

    return I

def circ_segment_void_area_avg(cond, param:str, hstar:float, int_err = 10**-4) -> float:
    """Void-weighted area-average over a circular segment defined by :math:`h^{*}`

    *WARNING*: uses ``scipy.integrate.dblquad``. May be computationally expensive
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        - ``hstar``: :math:`h^{*}` defining the circular segment of interest. See :any:`find_hstar_pos`
        - ``int_err``: acceptable warning for integral error (passed to scipy.integrate.dblquad). Defaults to 10**-4.
    
    **Raises**:
    
        - ``KeyError``: If ``param`` not available
    
    **Returns**:
    
        - Void-weighted area-average over the circular segment
    """

    # Check that the parameter that the user requested exists
    if not cond.check_param(param):
        raise KeyError(f"Invalid parameter {param} selected. Not present at {cond.check_param_loc(param)}")
    
    
    cond.mirror()

    rs = []
    phis = []
    vals = []
    denom = []
    for angle, r_dict in cond.data.items():
        for r_star, midas_dict in r_dict.items():
            rs.append(r_star)
            phis.append(angle * np.pi/180) # Convert degrees to radians
            denom.append(r_star * midas_dict['alpha'])
            try:
                vals.append(midas_dict[param] * r_star * midas_dict['alpha']) # Don't forget r dr dθ
            except:
                vals.append(np.NaN * r_star) # Don't forget r dr dθ
                print(f"Could not find {param} for φ = {angle}, r = {r_star}. Substituting NaN")

    # Set up interpolation
    rs = np.asarray(rs)
    phis = np.asarray(phis)
    vals = np.asarray(vals)   
    denom = np.asarray(denom)   

    interp = interpolate.LinearNDInterpolator(list(zip(rs, phis)), vals)
    interp2 = interpolate.LinearNDInterpolator(list(zip(rs, phis)), denom)

    # To integrate, need to get the bounds. They're different depending on if h* > 1
    def lower_r_bound_pos(phi):
        return max((1 - hstar) / np.sin(phi), 0)

    def upper_r_bound_neg(phi):
        if (phi <= 3*np.pi/2 - np.arccos(hstar-1)) or (phi >= 3*np.pi/2 + np.arccos(hstar-1)):
            return 1
        else:
            return (1 - hstar) / np.sin(phi)

    if hstar <= 1:
        I = integrate.dblquad(interp, np.pi/2 - np.arccos(1-hstar), np.pi/2 + np.arccos(1-hstar), lower_r_bound_pos, 1, epsabs=int_err )[0] / integrate.dblquad(interp2, np.pi/2 - np.arccos(1-hstar), np.pi/2 + np.arccos(1-hstar), lower_r_bound_pos, 1, epsabs=int_err )[0]
    elif hstar > 1:
        I = integrate.dblquad(interp, 0, 2*np.pi, 0, upper_r_bound_neg, epsabs=int_err )[0] / integrate.dblquad(interp2, 0, 2*np.pi, 0, upper_r_bound_neg, epsabs=int_err )[0]

    return I

def line_avg(cond, param:str, phi_angle:float) -> float:
    """Method for calculating the average value of param across a diameter defined by :math:`\\varphi=\\varphi^{*}` = angle

    .. math:: \\langle \\psi \\rangle_{L} = \\frac{1}{2} \int_L \\psi(r,\\varphi^{*}) \,dL
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        - ``phi_angle``: angle to calculate line average across
    
    **Raises**:
    
        - ``KeyError``: If ``param`` not available
    
    **Returns**:
    
        - Line average
    """

    # Check that the parameter that the user requested exists
    cond.mirror()

    if phi_angle not in cond.data.keys():
        if debug: print(cond.data, file=debugFID)
        print(f"Could not area-average {param} for condition {cond.name}\nData for {phi_angle} not found after mirroring!")
        return


    if not cond.check_param(param):
        raise KeyError(f"Invalid parameter {param} selected. Not present at {cond.check_param_loc(param)}")

    r_for_int = []
    var_for_int = []

    for rstar, midas_dict in cond.data[phi_angle].items():
        if rstar not in r_for_int:
            r_for_int.append(rstar)
            var_for_int.append(midas_dict[param])

    if phi_angle <=180:
        comp_angle = phi_angle+180
        
    else:
        comp_angle = phi_angle - 180
    
    for rstar, midas_dict in cond.data[comp_angle].items():
        if rstar not in r_for_int:
            r_for_int.append(-rstar)
            var_for_int.append(midas_dict[param])

    var_for_int = [param for _, param in sorted(zip(r_for_int, var_for_int))]
    r_for_int = sorted(r_for_int)

    I = integrate.simpson(y=var_for_int, x=r_for_int) / 2 # Integrate wrt theta, divide by normalized length

    return I

def line_avg_dev(cond, param:str, phi_angle:float, even_opt='first') -> float:
    """Method for calculating the average value of param across a diameter defined by :math:`\\varphi` = angle

    .. math:: \\langle \\psi \\rangle_{L} = \\frac{\int_L (\\psi(r,\\varphi^{*}) - \\langle \\psi \\rangle)^{2} \,dL}{\\langle \\psi \\rangle} 

    **Args**:
    
        - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        - ``phi_angle``: angle to calculate line average across
    
    **Raises**:
    
        - ``KeyError``: If ``param`` not available
    
    **Returns**:
    
        - Line average deviation
    """

    cond.mirror()

    if phi_angle not in cond.data.keys():
        if debug: print(cond.data, file=debugFID)
        print(f"Could not area-average {param} for condition {cond.name}\nData for {phi_angle} not found after mirroring!")
        return


    if not cond.check_param(param):
        raise KeyError(f"Invalid parameter {param} selected. Not present at {cond.check_param_loc(param)}")


    r_for_int = []
    var_for_int = []

    for rstar, midas_dict in cond.data[phi_angle].items():
        if rstar not in r_for_int:
            r_for_int.append(rstar)
            var_for_int.append((midas_dict[param] - area_avg(cond,param))**2)

    if phi_angle <=180:
        comp_angle = phi_angle+180
        
    else:
        comp_angle = phi_angle - 180
    
    for rstar, midas_dict in cond.data[comp_angle].items():
        if rstar not in r_for_int:
            r_for_int.append(-rstar)
            var_for_int.append((midas_dict[param] - area_avg(cond,param))**2)

    var_for_int = [param for _, param in sorted(zip(r_for_int, var_for_int))]
    r_for_int = sorted(r_for_int)

    I = integrate.simpson(y=var_for_int, x=r_for_int) / 2 / area_avg(cond,param)**2 # Integrate wrt theta, divide by normalized length

    return I

def void_area_avg(cond, param: str, even_opt='first', method = None) -> float:
    """Method for calculating the void-weighted area-average of a parameter

    .. math:: \\langle \\langle \\psi \\rangle \\rangle = \\frac{\iint_A \\alpha \\psi(r,\\varphi) r \,dr\,d\\varphi}{\iint_A \\alpha r \,dr\,d\\varphi} 
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter to void-weighted area-average. See :func:`~MARIGOLD.Condition.print_params` for options
        - ``even_opt``: option for ``scipy.integrate.simpsons``. Defaults to 'first'.
        - ``recalc``: if true, will recalculate area average. Defaults to True.
        - ``method``: method to void-weighted area-average. Defaults to None.

            - ``legacy``, using the same method as the Excel spreadsheets
            - ``legacy_old``, actually what we use in spreadsheets, hardcoded values
            - None, will use ``scipy.integrate.simpsons``. Recommended option
    
    **Raises**:
    
        - ``KeyError``: if ``param`` not found
    
    **Returns**:
    
        - void-weighted area-averaged value
    """

    # Check that the parameter that the user requested exists
    if not cond.check_param(param):
        raise KeyError(f"Invalid parameter {param} selected. Not present at {cond.check_param_loc(param)}")


    if method == 'legacy':
        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        if not cond.mirrored:
            warnings.warn("Mirroring in area-avg")
            cond.mirror()

        for angle, r_dict in cond.data.items():

            if angle == 360:    # We already have 0
                continue

            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi / 180)
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    try:
                        rs_temp.append(rstar)
                        vars_temp.append(float(midas_dict[param] * midas_dict['alpha'] * rstar))
                    except:
                        if debug: print('Problem with:', angle, r_dict, param)
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            vars.reverse()      # I noticed too late that the list goes from r/R = 0 to 1, not the other way around
            rs.reverse()        # Already wrote the actual Simpson's rule part, easier to just do this

            if debug: print("Arrays to integrate", angle, rs, vars, file=debugFID)

            if len(rs) != len(vars):
                ValueError( f"rs to integrate over {rs} must be the same length as params {vars}, occured at {angle}" )

            delta = abs(np.diff(rs))                # r/R steps

            la = 0
            for idx, var in enumerate(vars):
                coeff = 2 * (2**(idx % 2)) / 3      # Simpson's Rule coefficient
                
                if idx < len(delta):
                    if idx > 0:
                        if delta[idx - 1] == delta[idx]:
                            S = delta[idx] * coeff * var
                        else:
                            coeff = coeff / 2
                            S = (delta[idx - 1] * coeff * var) + (delta[idx] * coeff * var)
                    else:
                        S = delta[idx] * coeff * var
                else:
                    S = delta[idx - 1] * coeff * var
                
                la = la + S
            
            param_r.append(la)

        I = sum(param_r) / 8 / area_avg(cond,'alpha',method=method)

    elif method == 'legacy_old':
        # We have to integrate twice, once with resepect to r, again with respect to phi
        # Start with r

        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        if not cond.mirrored:
            warnings.warn("Mirroring in area-avg")
            cond.mirror()

        for angle, r_dict in cond.data.items():
            if angle == 360:
                continue

            try:
                if 0.95 in r_dict:
                    S_1 = 0.05 * sum((1 * 0,
                            4 * 0.95 * r_dict[0.95][param] * r_dict[0.95]['alpha'],
                            2 * 0.90 * r_dict[0.90][param] * r_dict[0.90]['alpha'],
                            4 * 0.85 * r_dict[0.85][param] * r_dict[0.85]['alpha'],
                            1 * 0.80 * r_dict[0.80][param] * r_dict[0.80]['alpha'],
                            )) / 3
                    
                    S_2 = 0.10 * sum((1 * 0.80 * r_dict[0.8][param] * r_dict[0.80]['alpha'],
                            4 * 0.70 * r_dict[0.7][param] * r_dict[0.7]['alpha'],
                            2 * 0.60 * r_dict[0.6][param] * r_dict[0.6]['alpha'],
                            4 * 0.50 * r_dict[0.5][param] * r_dict[0.5]['alpha'],
                            1 * 0.40 * r_dict[0.4][param] * r_dict[0.4]['alpha'],
                            )) / 3
                    
                    S_3 = 0.20 * sum((1 * 0.40 * r_dict[0.4][param] * r_dict[0.4]['alpha'],
                            4 * 0.20 * r_dict[0.2][param] * r_dict[0.2]['alpha'],
                            2 * 0.00 * r_dict[0.0][param] * r_dict[0.0]['alpha'],
                            )) / 3
                else:
                    S_1 = 0

                    S_2 = 0.10 * sum((1 * 0,
                            4 * 0.90 * r_dict[0.9][param] * r_dict[0.9]['alpha'],
                            2 * 0.80 * r_dict[0.8][param] * r_dict[0.8]['alpha'],
                            4 * 0.70 * r_dict[0.7][param] * r_dict[0.7]['alpha'],
                            2 * 0.60 * r_dict[0.6][param] * r_dict[0.6]['alpha'],
                            4 * 0.50 * r_dict[0.5][param] * r_dict[0.5]['alpha'],
                            1 * 0.40 * r_dict[0.4][param] * r_dict[0.4]['alpha'],
                            )) / 3
                    
                    S_3 = 0.20 * sum((1 * 0.40 * r_dict[0.4][param] * r_dict[0.4]['alpha'],
                            4 * 0.20 * r_dict[0.2][param] * r_dict[0.2]['alpha'],
                            2 * 0.00 * r_dict[0.0][param] * r_dict[0.0]['alpha'],               # Might be doubling up
                            )) / 3

                param_r.append(sum((S_1, S_2, S_3)))

            except Exception as e:
                print(e)
                
        I = sum(param_r) / 8 / area_avg(cond,'alpha',method=method)

    else:
        # We have to integrate twice, once with resepect to r, again with respect to phi
        # Start with r

        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        cond.mirror()


        for angle, r_dict in cond.data.items():

            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0: # This should be unnecessary now with the new mirror, but it's not hurting anyone by being here
                    try:
                        rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                        vars_temp.append( midas_dict[param] * midas_dict['alpha'] * rstar)
                    except:
                        if debug: print('Problem with:', angle, r_dict, param)
                    #if debug: print(angle, midas_dict, file=debugFID)
            
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            if debug: print("Arrays to integrate", angle, rs, vars, file=debugFID)
                
            param_r.append( integrate.simpson(y=vars, x=rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(y=vars, x=rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)

        param_r = [param for _, param in sorted(zip(angles, param_r))]
        angles = sorted(angles)

        I = integrate.simpson(y=param_r, x=angles) / np.pi / area_avg(cond,'alpha') # Integrate wrt theta, divide by normalized area

    return I

def interp_area_avg(cond, param:str, interp_type = 'linear', int_error = 10**-6) -> float:
    """Method for calculating the area-average of a parameter, based on an interpolation of the data
    
    *WARNING*: uses ``scipy.integrate.nquad``. May be computationally expensive
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        - ``interp_type``: interpolation type. See :any:`__call__`. Defaults to 'linear'.
        - ``int_error``: acceptable warning for integral error (passed to scipy.integrate.nquad). Defaults to 10**-6.
    
    **Returns**:
    
        - area-average of param
    """

    def integrand(phi, r): # phi will be in radians from dblquad
        return cond(phi, r, param, interp_method=interp_type) * r
    
    I = integrate.nquad(integrand, ranges = [(0, 1), (0, np.pi * 2)], opts = {'epsabs': int_error, 'points': list(cond.data[0].keys()), 'points': list(cond.data.keys())})[0] / cond.interp_area_avg('alpha') / np.pi
    return I

def interp_void_area_avg(cond, param:str, interp_type = 'linear', int_error = 10**-6) -> float:
    """Method for calculating the void-weighted area-average of a parameter, based on an interpolation of the data
    
    *WARNING*: uses ``scipy.integrate.nquad``. May be computationally expensive
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        - ``interp_type``: interpolation type. See :any:`__call__`. Defaults to 'linear'.
        - ``int_error``: acceptable warning for integral error (passed to scipy.integrate.nquad).. Defaults to 10**-6.
    
    **Returns**:
    
        - void-weighted area-average of param
    """

    def integrand(phi, r): # phi will be in radians from dblquad
        return cond(phi, r, param, interp_method=interp_type) * cond(phi, r, 'alpha', interp_method=interp_type) * r
    
    I = integrate.nquad(integrand, ranges = [(0, 1), (0, np.pi * 2)], opts = {'epsabs': int_error, 'points': list(cond.data.keys()), 'points': list(cond.data[0].keys())} )[0] / cond.interp_area_avg('alpha') / np.pi
    return I

def spline_void_area_avg(cond, param:str, int_error = 10**-6) -> float:
    """Function to void-weighted area-average param based on a spline interpolation
    
    *WARNING*: uses ``scipy.integrate.dblquad``. May be computationally expensive

    **Args**:
    
        - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        - ``int_error``: acceptable warning for integral error (passed to scipy.integrate.dblquad). Defaults to 10**-6.
    
    **Returns**:
    
        - integrand result
    """

    def integrand(phi, r):
        return cond.spline_interp[param](phi * 180/np.pi, r) * cond.spline_interp['alpha'](phi * 180/np.pi, r) * r
    
    def integrand_denom(phi, r):
        return cond.spline_interp['alpha'](phi * 180/np.pi, r) * r
    
    I = integrate.dblquad(integrand, 0, 1, 0, np.pi * 2, epsabs = int_error)[0] / integrate.dblquad(integrand_denom, 0, 1, 0, np.pi * 2, epsabs = int_error)[0]
    return I

def spline_circ_seg_area_avg(cond, param:str, hstar:float, int_err = 10**-4) -> float:
    """Function to area-average over a circular segment using the spline interpolation of param

    *WARNING*: uses ``scipy.integrate.dblquad``. May be computationally expensive
    
    **Args**:
    
        - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        - ``hstar``: distance from the top of the pipe that defines circular segment
        - ``int_err``: acceptable warning for integral error (passed to scipy.integrate.dblquad). Defaults to 10**-4.
    
    **Returns**:
    
        - integrand result
    """
    # Function to area-average over a circular segment using the spline interpolation of param

    # Inputs:
    #  - param, string of local parameter to area-average
    #  - hstar, distance from the top of the pipe that defines circular segment
    #  - int_err, acceptable warning for integral error (passed to scipy.integrate.dblquad)

    # Returns:
    #  - integrand result

    # Uses scipy.integrate.dblquad. May be computationally expensive
    
    # 

    def integrand(r, phi):
        return cond.spline_interp[param](phi * 180/np.pi, r) * r

    def integrand_denom(r, phi):
        return r

    def lower_r_bound_pos(phi):
        return max((1 - hstar) / np.sin(phi), 0)

    def upper_r_bound_neg(phi):
        if (phi <= 3*np.pi/2 - np.arccos(hstar-1)) or (phi >= 3*np.pi/2 + np.arccos(hstar-1)):
            return 1
        else:
            return (1 - hstar) / np.sin(phi)

    if hstar <= 1:
        I = integrate.dblquad(integrand, np.pi/2 - np.arccos(1-hstar), np.pi/2 + np.arccos(1-hstar), lower_r_bound_pos, 1, epsabs=int_err )[0] / integrate.dblquad(integrand_denom, np.pi/2 - np.arccos(1-hstar), np.pi/2 + np.arccos(1-hstar), lower_r_bound_pos, 1, epsabs=int_err )[0]
    elif hstar > 1:
        I = integrate.dblquad(integrand, 0, 2*np.pi, 0, upper_r_bound_neg, epsabs=int_err )[0] / integrate.dblquad(integrand_denom, 0, 2*np.pi, 0, upper_r_bound_neg, epsabs=int_err )[0]

    return I
