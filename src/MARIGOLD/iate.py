from .config import *
from .operations import *
from .correlations import *
import warnings

############################################################################################################################
#                                                                                                                          #
#                                                          IATE                                                            #
#                                                                                                                          #
############################################################################################################################
def iate_1d_1g(
        # Basic inputs
        cond, query, z_step = 0.01, io = None, geometry = None, R_c = None, L_res = None, cond2 = None,
        
        # IATE Coefficients
        C_WE = None, C_RC = None, C_TI = None, alpha_max = 0.75, C = 3, We_cr = 6, acrit_flag = 0, acrit = 0.13,

        # Method arguments
        preset = None, avg_method = None, cov_method = 'fixed', reconstruct_flag = False, cd_method = 'doe', dpdz_method = 'LM', void_method = 'driftflux',

        # Covariance calculation
        COV_RC = None, COV_TI = None,

        # Pressure drop calculation
        LM_C = 25, k_m = 0.10, m = 0.316, n = 0.25,

        # Void fraction calculation
        C0 = None, C_inf = 1.20,

        # Debugging
        verbose = False,

        ):
    """_summary_
    
    **Args**:
    
     - ``cond``: Condition object, part of MARIGOLD framework
     - ``query``: L/D endpoint
     - ``z_step``: Axial mesh cell size [-]. Defaults to 0.01.
     - ``R_c``: Radius of curvature ratio. Defaults to None.
     - ``io``: Output package of iate_1d_1g(), can be used as input for subsequent runs. Defaults to None.
     - ``geometry``: Geometry type, defaults to None, can be set to 'elbow', 'ubend'. Defaults to None.
     - ``cond2``: Second condition object, for possible interpolation. Defaults to None.
     - ``C_WE``: Wake entrainment coefficient. Defaults to None.
     - ``C_RC``: Random collision coefficient. Defaults to None.
     - ``C_TI``: Turbulent impact coefficient. Defaults to None.
     - ``alpha_max``: Maximum void fraction based on HCP bubble distribution, used for random collision calculation. Defaults to 0.75.
     - ``C``: Additional factor accounting for the range of eddy size capable of transporting bubbles. Value assumed to be 3, but no justification for this selection is made. Defaults to 3.
     - ``We_cr``: Weber number criterion, used for turbulent impact calculation. Defaults to 6.
     - ``acrit_flag``: Enable/disable shutting off turbulence-based mechanisms beyond a critical void fraction. Defaults to 0.
     - ``acrit``: Critical void fraction for shutting off turbulence-based mechanisms. Defaults to 0.13.
     - ``preset``: Author preset, fixes coefficients and method arguments to match old MATLAB runs. Defaults to None.
     - ``avg_method``: Area-averaging method, can be set to None (for Python Simpson's rule), 'legacy' (for Excel Simpson's Rule). Defaults to None.
     - ``cov_method``: Covariance calculation method. Defaults to 'fixed'.
     - ``reconstruct_flag``: Void reconstruction flag. Defaults to False.
     - ``cd_method``: Drag coefficient prediction method, 'err_iter', 'fixed_iter', or 'doe'. Defaults to 'doe'.
     - ``dpdz_method``: Pressure drop prediction method, 'LM' or 'Kim'. Defaults to 'LM'.
     - ``void_method``: Void fraction prediction method, 'driftflux' or 'continuity'. Defaults to 'driftflux'.
     - ``LM_C``: Lockhart-Martinelli Chisholm parameter. Defaults to 25.
     - ``k_m``: Minor loss coefficient. Defaults to 0.10.
     - ``L_res``: Restriction length. Defaults to None.
     - ``COV_RC``: Random collision mechanism covariance. Defaults to None.
     - ``COV_TI``: Turbulent impact mechanism covariance. Defaults to None.
     - ``m``: Friction factor constant. Defaults to 0.316.
     - ``n``: Friction factor constant. Defaults to 0.25.
     - ``C_inf``: Drift flux distribution parameter limiting value. Will be used to calculate C0, if none specified. Defaults to 1.20.
     - ``verbose``: Print warnings and commentary. Defaults to False.
    
    **Raises**:
    
     - ``ValueError``: _description_
    """
    # Code Repair:
    #  - Investigate L_res uses. May be replaceable with mesh length. Must be something Zhengting added.
    #  - There is no reason for jg_loc to be an array. It was fine before, don't know why Zhengting did this.
    #  - Couldn't handle it anymore. Removed Zhengting's edits; will need to go back and re-implement U-bend dissipation modeling approach.

    # Notes:
    #  - Notice some grav terms are made absolute; needs downward flow fixes
    #  - IATE coefficients set as optional inputs, with default values set depending on geometry
    #  - vgz calculation in elbow and dissipation length regions still need to be implemented
    #  - Need a way to compute void fraction across restrictions, void fraction prediction falters
    #  - Modify MG for Yadav data extraction
    #  - Revise vgj calculation

    # Apply IATE preset values
    preset_args = apply_preset(preset, cond)
    locals().update(preset_args)

    # MARIGOLD retrieval and setup
    theta           = cond.theta                                # Pipe inclination angle [degrees]
    Dh              = cond.Dh                                   # Hydraulic diameter [m]
    rho_f           = cond.rho_f                                # Liquid phase density [kg/m**3]
    rho_g           = cond.rho_g                                # Gas phase density [kg/m**3] or cond.rho_g 
    mu_f            = cond.mu_f                                 # Dynamic viscosity of water [Pa-s]
    mu_g            = cond.mu_g                                 # Dynamic viscosity of air [Pa-s]
    sigma           = cond.sigma                                # Surface tension of air/water [N/m]
    grav            = cond.gz                                   # Gravity constant (added by Drew to account for pipe inclination) (also, notably negative for VD)
    p_atm           = cond.p_atm                                # Ambient pressure 101325 [Pa]
    R_spec          = cond.R_spec                               # Specific gas constant for dry air [J/kg-K]
    T               = cond.T                                    # Ambient absolute temperature [K], for calculating air density as a function of pressure along channel

    # Starting L/D
    if io == None:
        LoverD      = cond.LoverD                               # Condition L/D
    else:
        LoverD      = io["z_mesh"][-1] / Dh                     # Last L/D
    
    # Mesh generation
    if query < LoverD+z_step:
        raise ValueError('Please choose a query L/D downstream of the boundary condition.')
    
    z_mesh = np.arange(LoverD, query + z_step, z_step)          # Axial mesh [-]
    z_mesh = z_mesh * Dh                                        # Axial mesh [m], units necessary for dp calculation
    z_step = z_step * Dh

    if geometry != None:
        R_c = R_c * Dh                                              # Radius of curvature

    ############################################################################################################################
    #                                                                                                                          #
    #                                                      COEFFICIENTS                                                        #
    #                                                                                                                          #
    ############################################################################################################################

    default_C_WE, default_C_RC, default_C_TI, overrides = set_geometry_coeffs(theta, geometry, Dh, R_c)
    locals().update(overrides)

    C_WE = C_WE if C_WE is not None else default_C_WE
    C_RC = C_RC if C_RC is not None else default_C_RC
    C_TI = C_TI if C_TI is not None else default_C_TI

    ############################################################################################################################
    #                                                                                                                          #
    #                                                       COVARIANCE                                                         #
    #                                                                                                                          #
    ############################################################################################################################
    if cond2 == None:
        cov_method = 'fixed'

    if cov_method == 'interp':
        # Use data at initial condition, void reconstruction downstream
        rf1, rf2 = (False, reconstruct_flag) if io is None else (reconstruct_flag, reconstruct_flag)
        
        if COV_RC == None:
            COV_RC1 = np.nan_to_num(cond.calc_COV_RC(reconstruct_flag = rf1, avg_method = avg_method, debug = False), nan=1.0)
            COV_RC2 = np.nan_to_num(cond2.calc_COV_RC(reconstruct_flag = rf2, avg_method = avg_method, debug = False), nan=1.0)
            COV_RC = np.interp(z_mesh / Dh,(cond.LoverD, cond2.LoverD),(COV_RC1, COV_RC2))
        
        if COV_TI == None:
            COV_TI1 = np.nan_to_num(cond.calc_COV_TI(reconstruct_flag = rf1, avg_method = avg_method, We_cr = We_cr, debug = False), nan=1.0)
            COV_TI2 = np.nan_to_num(cond2.calc_COV_TI(reconstruct_flag = rf2, avg_method = avg_method, We_cr = We_cr, debug = False), nan=1.0)
            COV_TI = np.interp(z_mesh / Dh,(cond.LoverD, cond2.LoverD),(COV_TI1, COV_TI2))

    else:
        if COV_RC == None:
            COV_RC = [1 for _ in range(len(z_mesh))]
        else:
            COV_RC = COV_RC * [1 for _ in range(len(z_mesh))]       # Not sure if this is necessary

        if COV_TI == None:
            COV_TI = [1 for _ in range(len(z_mesh))]
        else:
            COV_TI = COV_TI * [1 for _ in range(len(z_mesh))]

    ############################################################################################################################
    #                                                                                                                          #
    #                                                   BOUNDARY CONDITIONS                                                    #
    #                                                                                                                          #
    ############################################################################################################################
    
    # Variable initialization
    ai              = np.empty(len(z_mesh))
    alpha           = np.empty(len(z_mesh))
    Db              = np.empty(len(z_mesh))
    vgz             = np.empty(len(z_mesh))
    
    SWE             = np.empty(len(z_mesh))
    SRC             = np.empty(len(z_mesh))
    STI             = np.empty(len(z_mesh))
    SEXP            = np.empty(len(z_mesh))
    SVG             = np.empty(len(z_mesh))
    
    aiwe            = np.empty(len(z_mesh))
    airc            = np.empty(len(z_mesh))
    aiti            = np.empty(len(z_mesh))
    aiexp           = np.empty(len(z_mesh))
    aivg            = np.empty(len(z_mesh))

    if io == None:
        aiwe[0]     = 0
        airc[0]     = 0
        aiti[0]     = 0
        aiexp[0]    = 0
        aivg[0]     = 0

        jf          = cond.jf                                   # [m/s]
        jgloc       = cond.jgloc                                # [m/s]
        jgatm       = cond.jgatm                                # [m/s]

        if preset == 'kim':
            # jgloc = cond.jgref      # Testing for Bettis data

            ai[0]       = cond.area_avg_ai_sheet
            alpha[0]    = cond.area_avg_void_sheet
            Db[0]       = 6 * alpha[0] / ai[0]
        else:
            ai[0]       = area_avg(cond,"ai",method=avg_method)                     # [1/m]
            alpha[0]    = area_avg(cond,"alpha",method=avg_method)                  # [-]
            Db[0]       = void_area_avg(cond,"Dsm1",method=avg_method) / 1000       # [m]

    else:
        aiwe[0]     = io["aiwe"][-1]
        airc[0]     = io["airc"][-1]
        aiti[0]     = io["aiti"][-1]
        aiexp[0]    = io["aiexp"][-1]
        aivg[0]     = io["aivg"][-1]

        ai[0]       = io["ai"][-1]
        alpha[0]    = io["alpha"][-1]
        Db[0]       = io["Db"][-1]
        jf          = io["jf"]
        jgloc       = io["jgloc"]
        jgatm       = io["jgatm"]

    ########################################################################################################################
    # Pressure drop [Pa/m]

    # Calculate height change for gravitational loss
    if geometry == 'elbow':
        delta_h = (z_mesh[-1] - z_mesh[0]) * 2 / np.pi          # The height of an elbow is going to be its radius

    elif geometry == 'ubend':
        delta_h = 0

    else:
        delta_h = (z_mesh[-1] - z_mesh[0])                      # Dissipation region is going to be the same as standard VU

    # Calculate initial pressure and pressure gradient
    p = jgatm * p_atm / jgloc                                   # Back-calculate local corrected absolute pressure

    if preset == 'kim':
        p = cond.pz                                             # Override
        dpdz = cond.dpdz

    elif dpdz_method == 'interp':
        dpdz = ((cond2.jgatm * p_atm / cond2.jgloc) - p) / (cond2.LoverD - LoverD)

    else:
        dpdz = calc_dpdz(
            cond, 
            method = dpdz_method, 
            chisholm = LM_C, 
            m = m, 
            n = n, 
            k_m = k_m, 
            L = (query - LoverD) * Dh
            ) + ((rho_f * grav * delta_h) / (z_mesh[-1] - z_mesh[0]))   # Pressure gradient from gravity

    pz = p * (1 - (z_mesh - z_mesh[0]) * (dpdz / p))
    
	# Local gas density along the test section
    if preset == 'worosz':
        rho_gz = pz / R_spec / T                                # Worosz, Ideal Gas Law
    else:
        rho_gz = rho_g * pz / p                                 # Talley
    
    ############################################################################################################################
    #                                                                                                                          #
    #                                                           IATE                                                           #
    #                                                                                                                          #
    ############################################################################################################################
    
    # Calculate ai(z) to evaluate the steady state one-dim one-group model
    for i, z in enumerate(z_mesh):
        if (i+1) >= len(z_mesh):
            break
            
        jgloc = jgatm * p_atm / pz[i]                           # Talley used jgP1 and pressure at P1 instead of jgatm and p_atm
        
        if void_method == 'vgz_talley':
            vgz[i] = 1.05 * (jf + jgloc) - 1.23                 # Talley 2012, Eq. 3-31
            alpha[i] = jgloc / vgz[i]
            Db[i] = 6 * alpha[i] / ai[i]

        if void_method == 'vgz_interp' and cond2 != None:
            vgz[i] = np.interp(z_mesh[i] / Dh,
                               (cond.LoverD, cond2.LoverD),
                               (void_area_avg(cond,'ug1',method=avg_method), cond2.void_area_avg('ug1',method=avg_method))
                               )

            alpha[i] = jgloc / vgz[i]
            Db[i] = 6 * alpha[i] / ai[i]
            
        else:
            vgz[i] = jgloc / alpha[i]                           # Estimate void weighted velocity

        vfz = jf / (1 - alpha[i])

        ########################################################################################################################
        # Estimate bubble relative velocity <ur> (See Talley, 2012, 4.2.2.6)
        ur = 0.2                                                # Set to constant value of 0.23 m/s by Schilling (2007)

        if cd_method == 'err_iter':
            err = 0.1
            while abs(err) > 0.000001:
                ReD = rho_f * ur * Db[i] * (1 - alpha[i]) / mu_f    # Bubble Reynolds number
                CDe = 24 * (1 + 0.1 * ReD**0.75) / ReD              # Drag coefficient

                # Relative velocity (Ishii and Chawla, 1979), implemented by Worosz accounting for density difference
                ure1 = (4 * abs(grav) * Db[i] * (rho_f - rho_gz[i]) * (1 - alpha[i]) / 3 / CDe / rho_f)**0.5
                err = (ure1 - ur) / ur

                ur = ure1
            ReD = rho_f * ur * Db[i] * (1 - alpha[i]) / mu_f        # Update bubble Reynolds number
            CDwe = 24 * (1 + 0.1 * ReD**0.75) / ReD                 # Update drag coefficient

        elif cd_method == 'doe':
            # Original DOE_MATLAB_IAC
            ReD = rho_f * ur * Db[i] * (1 - alpha[i]) / mu_f
            CDwe = 24 * (1 + 0.1 * ReD**0.75) / ReD
            ur = (4 * abs(grav) * Db[i] / 3 / CDwe)**0.5            # Interestingly, Yadav keeps 9.8 instead of changing grav for angle

        elif cd_method == 'fixed_iter':
            for loop_idx in range(25):
                ReD = rho_f * ur * Db[i] * (1 - alpha[i]) / mu_f
                CDwe = 24 * (1 + 0.1 * ReD**0.75) / ReD

                if preset == 'kim':
                    ur = (grav * Db[i] / 3 / CDwe)**0.5
                else:
                    ur = (4 * grav * Db[i] / 3 / CDwe)**0.5
            
            ReD = rho_f * ur * Db[i] * (1 - alpha[i]) / mu_f
            CDwe = 24 * (1 + 0.1 * ReD**0.75) / ReD

        ########################################################################################################################
        # Estimate Energy Dissipation Rate and Turbulent Velocity (See Talley, 2012, 4.2.2.3)
        #   > One-group models written using turbulent fluctuation velocity, while models implemented in TRACE are written using
        #     dissipation rate
    
        mu_m = mu_f / (1 - alpha[i])                            # Mixture viscosity, given by Ishii and Chawla (Eq. 4-10 in Kim, 1999)
        rho_m = (1 - alpha[i]) * rho_f + alpha[i] * rho_gz[i]   # Mixture density

        vm = (rho_f * (1 - alpha[i]) * vfz + rho_gz[i] * alpha[i] * vgz[i]) \
            / rho_m                                             # Mixture velocity
        
        if preset == 'kim':
            Rem = rho_f * vm * Dh / mu_m                        # Mixture Reynolds number
        else:
            Rem = rho_m * vm * Dh / mu_m                        # Yadav

        fTW = 0.316 * (mu_m / mu_f)**0.25 / Rem**0.25           # Two-phase friction factor
        e = fTW * (vm**3) / 2 / Dh                              # Energy dissipation rate (Wu et al., 1998; Kim, 1999)
        u_t = 1.4 * e**(1/3) * Db[i]**(1/3)                     # Turbulent velocity (Batchelor, 1951; Rotta, 1972)

        ########################################################################################################################
        # Estimate sources & sinks in the Interfacial Area Transport Eqn. (Part 1)

        # Sink due to Wake Entrainment
        if theta == 0 and geometry == None:
            SWE[i] = 0
        
        elif geometry == 'elbow':   # (WIP)
            SWE[i] = 0
            
        else:
            SWE[i] = C_WE * CDwe**(1/3) * ur * ai[i]**2 / 3 / np.pi
        
        # Sink due to Random Collisions
        RC1 = u_t * ai[i]**2 / alpha_max**(1/3) / (alpha_max**(1/3) - alpha[i]**(1/3))
        RC2 = 1 - np.exp(-C * alpha_max**(1/3) * alpha[i]**(1/3) / (alpha_max**(1/3) - alpha[i]**(1/3)))
        SRC[i] = COV_RC[i] * C_RC * RC1 * RC2 / 3 / np.pi
        
        # Source due to Turbulent Impact
        TI1 = u_t * ai[i]**2 / alpha[i]
        We = rho_f * u_t**2 * Db[i] / sigma                     # Weber number criterion

        if We > We_cr:
            TI2 = (1 - We_cr / We)**0.5 * np.exp(-We_cr / We)
        else:
            TI2 = 0
        STI[i] = COV_TI[i] * C_TI * TI1 * TI2 / 18

        ########################################################################################################################
        # Estimate sources & sinks in the Interfacial Area Transport Eqn. (Part 2)

        # Source due to Bubble Expansion
        if preset == 'kim' or preset == 'talley':
            SEXP[i] = -2 / 3 / pz[i] * ai[i] * vgz[i] * (-dpdz)     # Original DOE_MATLAB_IAC
        else:
            if i <= 2:      # Previously 3, but in MATLAB (1 indexing vs. 0 indexing)
                # Forward difference for first node
                SEXP[i] = -2 / 3 / rho_gz[i] * ai[i] * vgz[i] * (rho_gz[i+1] - rho_gz[i]) / z_step
            else:
                # Backwards difference for remaining nodes
                SEXP[i] = -2 / 3 / rho_gz[i] * ai[i] * vgz[i] * (rho_gz[i] - rho_gz[i-1]) / z_step

        # Source/sink due to Bubble Acceleration (advection in Yadav's script) (VG for velocity gradient)
        if i <= 2:
            dvg = 0
            dvgdz = 0

            # dvgdz = dvg / z_step                              # Gonna leave this here for a good chuckle
        else:
            dvg = vgz[i] - vgz[i-1]
            dvgdz = dvg / z_step

        SVG[i] = ai[i] * dvgdz

        if preset == 'kim':
            SVG[i] = 0                                          # Suppress VG for Bettis data

        # Critical void fraction to shut off turbulence-based mechanisms
        if acrit_flag == 1 and alpha[i] > acrit:
            STI[i] = 0
            SRC[i] = 0

        ########################################################################################################################
        # Estimate Interfacial Area Concentration at an axial location z
        ai[i+1]         = ai[i] + z_step * (STI[i] - SWE[i] - SRC[i] + SEXP[i] - SVG[i]) / vgz[i]
        aiti[i+1]       = aiti[i] + z_step * STI[i] / vgz[i]
        airc[i+1]       = airc[i] + z_step * SRC[i] / vgz[i]
        aiexp[i+1]      = aiexp[i] + z_step * SEXP[i] / vgz[i]
        aiwe[i+1]       = aiwe[i] + z_step * SWE[i] / vgz[i]
        aivg[i+1]       = aivg[i] + z_step * SVG[i] / vgz[i]

        ########################################################################################################################
        # Estimate Void Fraction for the next step calculation
        if void_method == 'driftflux':      # Drift Flux Model

            j = jgloc + jf
            
            # Drift Velocity
            # Applicable for void fractions less than 20%; for void fractions greater than 30%, use Kataoka and Ishii 1987 for drift-velocity
            vgj = (2**0.5) * (sigma * abs(grav) * (rho_f - rho_gz[i]) / (rho_f**2))**0.25 * (1 - alpha[i])**(1.75)
            
            if C0 == None:
                C0 = C_inf - (C_inf - 1) * np.sqrt(rho_gz[i]/rho_f)     # Round tube drift flux distribution parameter
            
            alpha[i+1] = jgloc / (C0 * j + vgj)

        elif void_method == 'continuity':   # Continuity
            # Original continuity method
            # alpha[i+1] = alpha[i] - alpha[i] / pz[i] * -dpdz * z_step

            # Yadav
            # No discernible difference between this method and original continuity method
            if i <= 2:
                # Specific form of continuity not involving velocity gradients to avoid starting issue
                alpha[i+1] = alpha[i] - alpha[i] * (rho_gz[i+1] - rho_gz[i]) / rho_gz[i]

            else:
                alpha[i+1] = alpha[i] - (alpha[i] / (rho_gz[i] * vgz[i])) * ((rho_gz[i] * vgz[i]) - (rho_gz[i-1] * vgz[i-1]))
        
        elif void_method == 'vgz_talley':
            # Alpha calculated in front, still going to double-calculate the i+1 step out of paranoia
            # The way Talley had this imnplemented in Model_Horz.m is weird

            jgloc = jgatm * p_atm / pz[i+1]

            vgz[i+1] = 1.05 * (jf + jgloc) - 1.23                 # Talley 2012, Eq. 3-31
            alpha[i+1] = jgloc / vgz[i+1]
            
        elif void_method == 'vgz_interp':
            # Alpha calculated in front, still going to double-calculate the i+1 step out of paranoia
            # The way Talley had this imnplemented in Model_Horz.m is weird

            jgloc = jgatm * p_atm / pz[i+1]

            vgz[i+1] = np.interp(z_mesh[i] / Dh,
                               (cond.LoverD, cond2.LoverD),
                               (void_area_avg(cond,'ug1',method=avg_method), cond2.void_area_avg('ug1',method=avg_method))
                               )
            alpha[i+1] = jgloc / vgz[i+1]

        elif void_method == 'pressure_kim':
            f_f, f_g = calc_fric(cond, m = m, n = n)
            dpdz_f = f_f * 1/Dh * rho_f * jf**2 / 2

            phi_f2 = (dpdz - ((rho_f * grav * delta_h) / (z_mesh[-1] - z_mesh[0]))) / dpdz_f

            rho_x = rho_gz[i] / rho_f
            mu_x = mu_g / mu_f
            L_x = query - LoverD             # Geometry length scale, = L/D_restriction
            Re_f = rho_f * jf * Dh / mu_f

            chiM_inv = (3.165 * k_m / L_x * Re_f**0.25)**0.5

            quad_A = 1
            quad_B = LM_C * (1 + (chiM_inv**2))**0.5
            quad_C = 1 + (chiM_inv**2) - phi_f2

            chi_inv = ((-quad_B + (quad_B**2 - 4 * quad_A * quad_C)**0.5) / (2 * quad_A))     # Quadratic formula to solve for 1/X
            alpha_x = (chi_inv**8 / rho_x**3 / mu_x)**(1/7)                                   # Solve for alpha/(1-alpha)

            alpha[i+1] = alpha_x / (alpha_x + 1)

        elif void_method == 'pressure_LM':
            f_f, f_g = calc_fric(cond, m = m, n = n)
            dpdz_f = f_f * 1/Dh * rho_f * jf**2 / 2

            phi_f2 = (dpdz - ((rho_f * grav * delta_h) / (z_mesh[-1] - z_mesh[0]))) / dpdz_f

            rho_x = rho_gz[i] / rho_f
            mu_x = mu_g / mu_f

            quad_A = 1
            quad_B = LM_C
            quad_C = 1 - phi_f2

            chi_inv = ((-quad_B + (quad_B**2 - 4 * quad_A * quad_C)**0.5) / (2 * quad_A))     # Quadratic formula to solve for 1/X
            alpha_x = (chi_inv**8 / rho_x**3 / mu_x)**(1/7)                                   # Solve for alpha/(1-alpha)

            alpha[i+1] = alpha_x / (alpha_x + 1)

        # Estimate Sauter mean diameter for the next step calculation
        Db[i+1] = 6 * alpha[i+1] / ai[i+1]
        
    io = {
        "ai"            : ai,
        "alpha"         : alpha,
        "Db"            : Db,
        "jf"            : jf,
        "jgloc"         : jgloc,
        "jgatm"         : jgatm,
        
        "aiti"          : aiti,
        "airc"          : airc,
        "aiexp"         : aiexp,
        "aiwe"          : aiwe,
        "aivg"          : aivg,
        "z_mesh"        : z_mesh,
        "pz"            : pz
    }

    return io

############################################################################################################################
#                                                                                                                          #
#                                                        MODULES                                                           #
#                                                                                                                          #
############################################################################################################################
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

############################################################################################################################
#                                                                                                                          #
#                                                       FORMULAS                                                           #
#                                                                                                                          #
############################################################################################################################
def lineq(x, m, x0, b):
    return float(max(m * (x - x0) + b, 0))

def quadeq(x, a, b, c):
    return float(max(), 0)

def calc_Re(rho, v, D, mu):
    return rho * v * D / mu

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