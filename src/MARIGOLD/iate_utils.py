from .config import *

def apply_preset(preset, cond):
    """Applies hard-coded presets to the condition object and returns overrides."""
    overrides = {}

    if preset == 'kim':
        cond.rho_f                  = 998
        cond.rho_g                  =  1.226
        cond.mu_f                   = 0.001
        cond.sigma                  = 0.07278
        cond.g                      = 9.8
        cond.gz                     = 9.8 * np.sin(np.radians(cond.theta))
        cond.p_atm                  = 101330

        overrides.update({
            "cd_method"             : 'fixed_iter',
            "C0"                    : 1.12
        })

    elif preset == 'talley':
        cond.rho_f                  = 998
        cond.rho_g                  = 1.23
        cond.mu_f                   = 0.001
        cond.mu_g                   = 1.73E-5
        cond.sigma                  = 0.07278
        cond.p_atm                  = 101353            # Equivalent to 14.7 [psi]

        overrides.update({
            "avg_method"            : 'legacy_old',
            "cd_method"             : 'doe',
            "dpdz_method"           : 'LM',
            "reconstruct_flag"      : True,
            "LM_C"                  : 25
        })

    elif preset == 'yadav':         # (WIP)
        overrides.update({
            "void_method"           : 'continuity'
        })

    elif preset == 'worosz':        # (WIP)
        overrides.update({
            "cd_method"             : 'err_iter'
        })

    elif preset == 'ryan':
        cond.Dh                     = 0.0254
        cond.p_atm                  = 101325
        cond.rho_f                  = 998
        cond.rho_g                  = 1.204
        cond.mu_f                   = 0.001002
        cond.sigma                  = 0.0728
        cond.R_spec                 = 287.058
        cond.T                      = 293.15

        # Double-check coefficients against thesis
        # interpolate.m opt = 2: linearly interpolate P, calculate <<vg>> from exp <jg>atm and P, calculate <alpha> using DF
        overrides.update({
            "z_step"                : 0.001,
            "alpha_max"             : 0.75,
            "C"                     : 3,
            "C_RC"                  : 0.004,
            "C_TI"                  : 0.085,
            "We_cr"                 : 6,
            "acrit_flag"            : 2,
            "acrit"                 : 0.13,
            "C_WE"                  : 0.004
        })

    elif preset == 'quan':
        cond.Dh                     = 0.0254
        cond.rho_f                  = 998
        cond.rho_g                  = 1.226
        cond.mu_f                   = 0.001
        cond.p_atm                  = 101353            # Equivalent to 14.7 [psi]

        overrides.update({
            "L_res"                 : 31.67,            # Length of restriction, based on U-bend experimental dpdz data
            "dpdz_method"           : 'kim',
            "C_WE"                  : 0.000,
            "C_RC"                  : 0.060,
            "C_TI"                  : 0.000,
            "acrit"                 : 1.00
        })

    return overrides

def set_geometry_coeffs(theta, geometry, Dh, R_c):
    """
    Assign IATE coefficients based on pipe geometry and orientation.

    For future reference, when adding other geometry types, want to maintain consistency in if/then logic.
     1. Check angles with straight pipes
     2. Check angles with restrictions
     3. Check restrictions
          a. 'elbow'
          b. 'ubend'
          c. 'dissipation'
          d. Other
     4. Else, default to vertical-upward
    
    Returns:
        - C_WE, C_RC, C_TI: IATE coefficients
        - overrides: additional parameter updates (e.g. We_cr, LM_C, k_m, acrit)
    """
    overrides = {}

    if theta == 0 and geometry is None:
        # Horizontal straight pipe (Talley, 2012)
        C_WE = 0.000
        C_RC = 0.003
        C_TI = 0.014

        overrides.update({
            "We_cr"                 : 5
        })

    elif geometry == 'elbow':
        # Elbow (Yadav, 2013)
        C_WE = 0.000
        C_RC = 0.008
        C_TI = 0.085

        overrides.update({
        })
        
    elif geometry == 'ubend':
        # U-bend (Quan, 2024)
        C_WE = 0.000
        C_RC = 0.010
        C_TI = 0.008

        overrides.update({
            "acrit"                 : 1.00,
            "We_cr"                 : 6,
            "LM_C"                  : 85,
            "k_m"                   : 0.20
        })

    elif geometry == 'dissipation':
        # U-bend dissipation region (Quan, 2025)
        C_WE = 0.000
        C_RC = 0.004
        C_TI = 0.085

        overrides.update({
            "acrit"                 : 1.00,
            "We_cr"                 : 6,
            "LM_C"                  : 68
        })

    elif geometry == 'vd':
        # Vertical-downward (Ishii et al.)
        C_WE = 0.002
        C_RC = 0.004
        C_TI = 0.034

        overrides.update({
        })

    else:
        # Default to vertical-upward (Ishii et al.)
        C_WE = 0.002
        C_RC = 0.004
        C_TI = 0.085

        overrides.update({
        })

    return C_WE, C_RC, C_TI, overrides

def print_iate_args(frame):
    """Print grouped arguments for iate_1d_1g based on current values."""
    values = inspect.getargvalues(frame).locals

    def print_group(title, arg_names):
        print(f"\n{title}:")
        for name in arg_names:
            print(f"  {name:<12} = {values[name]!r}")

    print_group("Basic inputs", ["cond", "query", "z_step", "io", "geometry", "R_c", "cond2"])
    print_group("IATE coefficients", ["C_WE", "C_RC", "C_TI", "alpha_max", "C", "We_cr", "acrit_flag", "acrit"])
    print_group("Method arguments", ["preset", "avg_method", "cov_method", "reconstruct_flag", "cd_method", "dpdz_method", "void_method"])
    print_group("Covariance calculation", ["COV_WE", "COV_RC", "COV_TI"])
    print_group("Pressure drop calculation", ["LM_C", "k_m", "m", "n"])
    print_group("Void fraction calculation", ["C0", "C_inf"])
    print_group("Debugging", ["verbose"])

############################################################################################################################
#                                                                                                                          #
#                                                      CORRELATIONS                                                        #
#                                                                                                                          #
############################################################################################################################
def calc_CD(Re, method = None):

    if method == '':
        CD = 24 * (1 + 0.1 * Re**0.75) / Re
        pass
    else:
        CD = 24 * (1 + 0.1 * Re**0.75) / Re
        pass

    return CD

def calc_ur():
    ur = 0
    return ur

def calc_void_dpdz(cond, jf, Dh, z_mesh, dpdz, LM_C, m, n, k_m, delta_h, grav, mu_f, mu_g, rho_f, rho_g, L_x, method):

    f_f, f_g = calc_fric(cond, m = m, n = n)
    dpdz_f = f_f * 1/Dh * rho_f * jf**2 / 2

    phi_f2 = (dpdz - ((rho_f * grav * delta_h) / (z_mesh[-1] - z_mesh[0]))) / dpdz_f

    rho_x = rho_g / rho_f
    mu_x = mu_g / mu_f

    if method == 'kim':
        Re_f = rho_f * jf * Dh / mu_f

        chiM_inv = (3.165 * k_m / L_x * Re_f**0.25)**0.5

        quad_A = 1
        quad_B = LM_C * (1 + (chiM_inv**2))**0.5
        quad_C = 1 + (chiM_inv**2) - phi_f2

    else:
        # Default to LM
        quad_A = 1
        quad_B = LM_C
        quad_C = 1 - phi_f2

    chi_inv = ((-quad_B + (quad_B**2 - 4 * quad_A * quad_C)**0.5) / (2 * quad_A))   # Quadratic formula to solve for 1/X
    alpha_x = (chi_inv**8 / rho_x**3 / mu_x)**(1/7)                                 # Solve for alpha/(1-alpha)

    alpha = alpha_x / (alpha_x + 1)

    return alpha

############################################################################################################################
#                                                                                                                          #
#                                                        METHODS                                                           #
#                                                                                                                          #
############################################################################################################################
def calc_COV(cond, alpha_peak = 0.75, alpha_cr = 0.11, We_cr = 5, avg_method = 'legacy', reconstruct_flag = True):
    """Calculates the experimental covariances based on Talley (2012) method (without modification factor m_RC)
        - Stored in cond.COV_XX
        
        Inputs:
        - alpha_peak, maximum void fraction based on hexagonal-closed-packed (HCP) bubble distribution
        - alpha_cr, critical alpha to activate Random Collision, Talley (2012), Kong (2018) 
    """

    if reconstruct_flag == True:
        alpha_str = 'alpha_reconstructed'
    else:
        alpha_str = 'alpha'

    # Temporary, replace later
    alpha_avg       = area_avg(cond,'alpha',method=avg_method)
    alpha_avg       = round(alpha_avg,3)

    ai_avg          = area_avg(cond,'ai',method=avg_method)
    ai_avg          = round(ai_avg,2)

    # Constants used by Talley were different in Excel from MATLAB
    rho_f           = cond.rho_f                                        # Liquid phase density [kg/m**3]
    rho_g           = cond.rho_g                                        # Gas phase density [kg/m**3]
    mu_f            = cond.mu_f                                         # Dynamic viscosity of water [Pa-s]
    sigma           = cond.sigma                                        # Surface tension
    Dh              = cond.Dh                                           # Hydraulic diameter [m]
            
    rho_m           = (1 - alpha_avg) * rho_f + alpha_avg * rho_g       # Mixture density
    mu_m            = mu_f / (1 - alpha_avg)                            # Mixture viscosity
    v_m             = (rho_f * cond.jf + rho_g * cond.jgloc) / rho_m    # Mixture velocity
    # Rem             = rho_m * v_m * Dh / mu_m                         # Mixture Reynolds number, ***CAREFUL*** I have seen some versions of the IATE script that use rho_f instead as an approximation
    Rem             = rho_m * v_m * Dh / mu_f                           # Mixture Reynolds number, older versions use mu_f as an approximation
    f_TP            = 0.316 * (mu_m / mu_f / Rem)**0.25                 # Two-phase friction factor, Talley (2012) and Worosz (2015), also used in iate_1d_1g
    eps             = f_TP * v_m**3 / 2 / Dh                            # Energy dissipation rate (Wu et al., 1998; Kim, 1999), also used in iate_1d_1g        
    eps             = round(eps,2)

    # Switch away from using data
    alpha_avg       = area_avg(cond,alpha_str,method=avg_method)

    Dsm_exp         = 1000 * 6 * alpha_avg / ai_avg
    Dsm_exp         = round(Dsm_exp,2) / 1000

    for angle, r_dict in cond.data.items():
        for rstar, midas_dict in r_dict.items():
            
            alpha_loc = midas_dict[alpha_str]
            
            if reconstruct_flag == True:
                Db_loc = Dsm_exp
                # ai_loc = ai_avg * alpha_loc / alpha_avg
                ai_loc = 6 * alpha_loc / Dsm_exp

                if ai_loc != 0:
                    Db_loc = 6 * alpha_loc / ai_loc
                else:
                    Db_loc = 0
            else:
                ai_loc = midas_dict['ai']
                # Db_loc = midas_dict['Dsm1'] / 1000
                
                if ai_loc != 0:
                    Db_loc = 6 * alpha_loc / ai_loc
                else:
                    Db_loc = 0
            
            if alpha_loc <= alpha_cr:                                   # Check if local void fraction is less than or equal to alpha_cr
                u_t = 1.4 * np.cbrt(eps) * np.cbrt(Db_loc)              # Turbulent velocity (Batchelor, 1951; Rotta, 1972), also used in iate_1d_1g
            else:
                u_t = 0                                                 # TI and RC are driven by the turbulent fluctuation velocity (u_t)

            # Talley 2012, section 3.3.1
            # COV_RC
            COV_RC_loc = u_t * ai_loc**2 / (np.cbrt(alpha_peak) * (np.cbrt(alpha_peak) - np.cbrt(alpha_loc)))
            midas_dict['COV_RC_loc'] = COV_RC_loc

            # COV_TI
            We = rho_f * u_t**2 * Db_loc / sigma                        # Weber number criterion
            if We >= We_cr:
                COV_TI_loc = (u_t * ai_loc**2 / alpha_loc) * np.sqrt(1 - (We_cr / We)) * np.exp(-We_cr / We)
            else:
                COV_TI_loc = 0
            midas_dict['COV_TI_loc'] = COV_TI_loc

            # Ryan 2022, section 7.2.2
            # COV_WE
            Reb = calc_Re()
            CD = calc_CD(Reb)
            ur = calc_ur()
            COV_WE_loc = np.cbrt(CD) * ai_loc**2 * ur
            midas_dict['COV_WE_loc'] = COV_WE_loc

    # Talley does not area-average local u_t; instead computes <u_t> with area-averaged parameters
    u_t_avg = 1.4 * np.cbrt(eps) * np.cbrt(6 * alpha_avg / ai_avg)
    We_avg = rho_f * u_t_avg**2 * (6 * alpha_avg / ai_avg) / sigma

    if u_t_avg > 0:
        COV_RC_avg = u_t_avg * ai_avg**2 / (np.cbrt(alpha_peak) * (np.cbrt(alpha_peak) - np.cbrt(alpha_avg)))
        COV_RC = area_avg(cond,'COV_RC_loc',method=avg_method) / COV_RC_avg

        COV_TI_avg = (u_t_avg * ai_avg**2 / alpha_avg) * np.sqrt(1 - (We_cr / We_avg)) * np.exp(-We_cr / We_avg)
        COV_TI = area_avg(cond,'COV_TI_loc',method=avg_method) / COV_TI_avg

    else:
        COV_RC = 0
        COV_TI = 0

    cond.COV_RC = COV_RC
    cond.COV_TI = COV_TI

    return COV_RC, COV_TI

def reconstruct_void(cond, method='talley', avg_method = 'legacy'):
    """Method to reconstruct the void fraction profile by various means. 
    
    **Args:**

        - ``method``: method to use to reconstruct void. Defaults to ``'talley'``. Options include:
            - ``'talley'``
            - ``'ryan'`` (WIP)
        - ``avg_method``: option passed to :func:`~MARIGOLD.Condition.Condition.area_avg`. Defaults to ``'legacy'``.
    
    Returns:
        - area-averaged reconstructed void
    """

    if method.lower() == 'talley':
        # STEP 1: Determine r* at which void fraction reaches zero.
        cond.roverRend = round(-1.472e-5 * cond.Ref + 2.571,1)      # Inner r/R, outer end fixed at r/R = 1. Also, Talley rounds his r/R_end to the nearest 0.1

        # STEP 3: Determine peak void fraction location. Peak location assumed to be 0.90.
        rstar_peak = 0.90

        # STEP 2: Determine peak void fraction value.
        def find_alpha_peak(alpha_peak, rstar_peak=0.90):
            # Talley's reconstruction has a finer grid than experimental data
            # This will affect the slope of the drop-off point if r/R_end is not coincident with a point on the experiment mesh
            r_points = np.arange(0,1,0.05)
            cond.add_mesh_points(r_points)

            interps = {}

            for angle, r_dict in cond.data.items():                 # 360 degrees covered, not just one quadrant
                
                if angle < 90:
                    interps[angle] = {'angle_nn': angle_q1, 'm_nn': m_i, 'x0_nn': rstar_anchor, 'b_nn': anchor}

                # Operate in first quadrant
                if angle <= 90:
                    angle_q1 = angle
                elif angle <= 180:
                    angle_q1 = 180 - angle
                elif angle <= 270:
                    angle_q1 = angle - 180
                else:
                    angle_q1 = 360 - angle

                # Save previous angle linear interpolations for nearest neighbor determination
                if angle_q1 != 90:
                    angle_nn = interps[angle_q1]['angle_nn']
                    m_nn = interps[angle_q1]['m_nn']
                    x0_nn = interps[angle_q1]['x0_nn']
                    b_nn = interps[angle_q1]['b_nn']

                # Find the peak nearest neighbor
                if angle_q1 == 90:
                    # For 90 degrees, just alpha_peak
                    peak = alpha_peak

                    anchor = 0
                    rstar_anchor = cond.roverRend

                elif angle_q1 == 0:
                    # For 0 degrees, value at r/R = 0 along the 90 degree axis
                    peak = lineq(x = 0, 
                                    m = 0 - alpha_peak / (cond.roverRend - rstar_peak),
                                    x0 = rstar_peak,
                                    b = alpha_peak)

                    anchor = peak
                    rstar_anchor = 0

                else:
                    # Find r* of previous angle at equivalent y-coordinate
                    y_peak = rstar_peak * np.sin(angle_q1 * np.pi / 180)
                    rstar_nn = y_peak / np.sin(angle_nn * np.pi / 180)

                    # Peak void fraction of current angle defined as void fraction at previous angle r*
                    peak = lineq(x = rstar_nn,
                                    m = m_nn,
                                    x0 = x0_nn,
                                    b = b_nn)
                    
                    anchor = 0
                    rstar_anchor = cond.roverRend / np.cos((90 - angle_q1) * np.pi / 180)       # Talley's implementation in Excel
                    # rstar_anchor = cond.roverRend / np.sin(angle_q1 * np.pi / 180)            # But why not like this
                    
                    # if cond.roverRend > 0:
                    #     anchor = 0
                    #     rstar_anchor = cond.roverRend
                    # else:
                    #     # Value at r/R = 0 along the 90 degree axis
                    #     anchor = lineq(x = 0,
                    #                    m = 0 - alpha_peak / (cond.roverRend - rstar_peak),
                    #                    x0 = rstar_peak,
                    #                    b = alpha_peak)
                    #     rstar_anchor = 0

                # Linear interpolation, y = mx + b
                m_o1 = (0 - peak) / (1 - rstar_peak)                    # Slope of outer interpolation
                m_i = (peak - anchor) / (rstar_peak - rstar_anchor)     # Slope of inner interpolation

                base = lineq(x = -rstar_peak,
                                m = m_i,
                                x0 = rstar_anchor,
                                b = anchor)
                rstar_base = -rstar_peak

                m_o2 = (0 - base) / (-1 - rstar_base)                   # Slope of outer interpolation, on opposite wall

                for rstar, midas_dict in r_dict.items():                # Only goes from 0.0 to 1.0
                    
                    if angle >= 180 and angle < 360:
                        rstar = -rstar                                  # Use same slopes for -1.0 to 0.0

                    # Alpha interpolations in terms of r* (see page 134 of Talley 2012)
                    if rstar > rstar_peak:
                        # Outer interpolation
                        midas_dict['alpha_reconstructed'] = lineq(x = rstar,
                                                                    m = m_o1,
                                                                    x0 = rstar_peak,
                                                                    b = peak)

                    elif abs(rstar) <= rstar_peak:
                        if angle_q1 == 0:
                            # Between peak locations, uniform profile
                            midas_dict['alpha_reconstructed'] = peak
                        else:
                            # Inner interpolation
                            midas_dict['alpha_reconstructed'] = lineq(x = rstar,
                                                                        m = m_i,
                                                                        x0 = rstar_anchor,
                                                                        b = anchor)
                            
                    else:
                        if angle_q1 == 0:
                            # Symmetry for 180 degrees
                            midas_dict['alpha_reconstructed'] = lineq(x = rstar,
                                                                        m = -m_o1,
                                                                        x0 = rstar_peak,
                                                                        b = peak)
                        else:
                            # Opposite wall interpolation
                            midas_dict['alpha_reconstructed'] = lineq(x = rstar,
                                                                        m = m_o2,
                                                                        x0 = rstar_base,
                                                                        b = base)

                cond.alpha_peak = alpha_peak
                
            return abs( round(area_avg(cond,'alpha',method=avg_method),3) - area_avg(cond,'alpha_reconstructed',method=avg_method) )    # Talley's value is to three decimals of precision

    elif method.lower() == 'ryan':
        # STEP 1: Determine r* at which void fraction reaches zero (Eq. 7-32).
        cond.roverRend = (1.3 - ((1.57e-5) * cond.Ref)) * np.cos(cond.theta * np.pi / 180)

        # STEP 2: Determine peak void fraction location (Eq. 7-31).
        if cond.Ref > 75000:
            rstar_peak = (1 - (2 * area_avg(cond,'Dsm') / cond.Dh)) * np.cos(cond.theta * np.pi / 180)**0.25
        else:
            rstar_peak = (1 - (2 * area_avg(cond,'Dsm') / cond.Dh))

        # STEP 3: Determine peak void fraction value (WIP).
        # This approach is nearly identical to Talley's, except the linear estimation lineq() is replaced with a quadratic quadeq().
        def find_alpha_peak(alpha_peak, rstar_peak=0.90):
            # Talley's reconstruction has a finer grid than experimental data
            # This will affect the slope of the drop-off point if r/R_end is not coincident with a point on the experiment mesh
            r_points = np.arange(0,1,0.05)
            cond.add_mesh_points(r_points)

            interps = {}
            
            for angle, r_dict in cond.data.items():                 # 360 degrees covered, not just one quadrant
                
                if angle < 90:
                    interps[angle] = {'angle_nn': angle_q1, 'm_nn': m_i, 'x0_nn': rstar_anchor, 'b_nn': anchor}

                # Operate in first quadrant
                if angle <= 90:
                    angle_q1 = angle
                elif angle <= 180:
                    angle_q1 = 180 - angle
                elif angle <= 270:
                    angle_q1 = angle - 180
                else:
                    angle_q1 = 360 - angle

                # Save previous angle linear interpolations for nearest neighbor determination
                if angle_q1 != 90:
                    angle_nn = interps[angle_q1]['angle_nn']
                    m_nn = interps[angle_q1]['m_nn']
                    x0_nn = interps[angle_q1]['x0_nn']
                    b_nn = interps[angle_q1]['b_nn']

                # Find the peak nearest neighbor
                if angle_q1 == 90:
                    # For 90 degrees, just alpha_peak
                    peak = alpha_peak

                    anchor = 0
                    rstar_anchor = cond.roverRend

                elif angle_q1 == 0:
                    # For 0 degrees, value at r/R = 0 along the 90 degree axis
                    peak = lineq(x = 0, 
                                    m = 0 - alpha_peak / (cond.roverRend - rstar_peak),
                                    x0 = rstar_peak,
                                    b = alpha_peak)

                    anchor = peak
                    rstar_anchor = 0

                else:
                    # Find r* of previous angle at equivalent y-coordinate
                    y_peak = rstar_peak * np.sin(angle_q1 * np.pi / 180)
                    rstar_nn = y_peak / np.sin(angle_nn * np.pi / 180)

                    # Peak void fraction of current angle defined as void fraction at previous angle r*
                    peak = lineq(x = rstar_nn,
                                    m = m_nn,
                                    x0 = x0_nn,
                                    b = b_nn)
                    
                    anchor = 0
                    rstar_anchor = cond.roverRend / np.cos((90 - angle_q1) * np.pi / 180)       # Talley's implementation in Excel
                    # rstar_anchor = cond.roverRend / np.sin(angle_q1 * np.pi / 180)            # But why not like this

                # Linear interpolation, y = mx + b
                m_o1 = (0 - peak) / (1 - rstar_peak)                    # Slope of outer interpolation
                m_i = (peak - anchor) / (rstar_peak - rstar_anchor)     # Slope of inner interpolation

                base = lineq(x = -rstar_peak,
                                m = m_i,
                                x0 = rstar_anchor,
                                b = anchor)
                rstar_base = -rstar_peak

                m_o2 = (0 - base) / (-1 - rstar_base)                   # Slope of outer interpolation, on opposite wall

                for rstar, midas_dict in r_dict.items():                # Only goes from 0.0 to 1.0

                    if angle >= 180 and angle < 360:
                        rstar = -rstar                                  # Use same slopes for -1.0 to 0.0

                    # Alpha interpolations in terms of r* (see page 134 of Talley 2012)
                    if rstar > rstar_peak:
                        # Outer interpolation
                        midas_dict['alpha_reconstructed'] = lineq(x = rstar,
                                                                    m = m_o1,
                                                                    x0 = rstar_peak,
                                                                    b = peak)

                    elif abs(rstar) <= rstar_peak:
                        if angle_q1 == 0:
                            # Between peak locations, uniform profile
                            midas_dict['alpha_reconstructed'] = peak
                        else:
                            # Inner interpolation
                            midas_dict['alpha_reconstructed'] = lineq(x = rstar,
                                                                        m = m_i,
                                                                        x0 = rstar_anchor,
                                                                        b = anchor)
                            
                    else:
                        if angle_q1 == 0:
                            # Symmetry for 180 degrees
                            midas_dict['alpha_reconstructed'] = lineq(x = rstar,
                                                                        m = -m_o1,
                                                                        x0 = rstar_peak,
                                                                        b = peak)
                        else:
                            # Opposite wall interpolation
                            midas_dict['alpha_reconstructed'] = lineq(x = rstar,
                                                                        m = m_o2,
                                                                        x0 = rstar_base,
                                                                        b = base)

                cond.alpha_peak = alpha_peak

            return abs( round(area_avg(cond,'alpha',method=avg_method),3) - area_avg(cond,'alpha_reconstructed',method=avg_method) )    # Talley's value is to three decimals of precision
    
    result = minimize(find_alpha_peak, x0=0.5, bounds=((0,1),))

    if result.success:
        cond.alpha_peak_reconstructed = result.x
        find_alpha_peak(cond.alpha_peak_reconstructed)
    else:
        warnings.warn("Minimization did not return a successful result")
        print(result.message)
    
    print(f"\tr/R_end: {cond.roverRend}")
    print(f"\talpha_peak: {cond.alpha_peak}")
    print(f"\talpha_data_aavg: {round(area_avg(cond,'alpha',method=avg_method),3)}")
    print(f"\talpha_reconstructed_aavg: {area_avg(cond,'alpha_reconstructed',method=avg_method)}")
    
    return area_avg(cond,"alpha_reconstructed")