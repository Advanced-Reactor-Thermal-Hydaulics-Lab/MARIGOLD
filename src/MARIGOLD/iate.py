from .config import *
from .iate_utils import *
from .operations import *

############################################################################################################################
#                                                                                                                          #
#                                                          IATE                                                            #
#                                                                                                                          #
############################################################################################################################
def iate_1d_1g(
        # Basic inputs
        cond, query, z_step = 0.01, io = None, geometry = None, R_c = None, cond2 = None,
        
        # IATE coefficients
        C_WE = None, C_RC = None, C_TI = None, alpha_max = 0.75, C = 3, We_cr = 6, acrit_flag = 0, acrit = 0.13,

        # Method arguments
        preset = None, avg_method = None, cov_method = 'fixed', reconstruct_flag = False, cd_method = 'doe', dpdz_method = 'LM', void_method = 'driftflux',

        # Covariance calculation
        COV_WE = None, COV_RC = None, COV_TI = None,

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
     - ``COV_WE``: Wake entrainment mechanism covariance. Defaults to None.
     - ``COV_RC``: Random collision mechanism covariance. Defaults to None.
     - ``COV_TI``: Turbulent impact mechanism covariance. Defaults to None.
     - ``m``: Friction factor constant. Defaults to 0.316.
     - ``n``: Friction factor constant. Defaults to 0.25.
     - ``C_inf``: Drift flux distribution parameter limiting value. Will be used to calculate C0, if none specified. Defaults to 1.20.
     - ``verbose``: Print warnings and commentary. Defaults to False.
    
    **Raises**:
    
     - ``ValueError``: _description_
    """
    # Notes:
    #  - IATE coefficients set as optional inputs, with default values set depending on geometry
    #  - Void fraction calculation methods are divided into three categories:
    #     - (1) interpolation-based methods, located before the IATE loop
    #     - (2) Prior step independent methods, located at the top of the IATE loop
    #     - (3) Prior step dependent methods, located at the end of the IATE loop

    # To do:
    #  - Notice some grav terms are made absolute; need downward flow fixes
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

    print_iate_args(inspect.currentframe())

    ############################################################################################################################
    #                                                                                                                          #
    #                                                       COVARIANCE                                                         #
    #                                                                                                                          #
    ############################################################################################################################

    if cov_method == 'interp' and cond2 != None:
        # Use data at initial condition, void reconstruction downstream
        rf1, rf2 = (False, reconstruct_flag) if io is None else (reconstruct_flag, reconstruct_flag)
        
        COV_WE1, COV_RC1, COV_TI1 = np.nan_to_num(
            calc_COV(cond, reconstruct_flag = rf1, avg_method = avg_method), nan=1.0
        )
        COV_WE2, COV_RC2, COV_TI2 = np.nan_to_num(
            calc_COV(cond2, reconstruct_flag = rf2, avg_method = avg_method), nan=1.0
        )

        COV_WE = np.interp(z_mesh / Dh,(cond.LoverD, cond2.LoverD),(COV_WE1, COV_WE2))
        COV_RC = np.interp(z_mesh / Dh,(cond.LoverD, cond2.LoverD),(COV_RC1, COV_RC2))
        COV_TI = np.interp(z_mesh / Dh,(cond.LoverD, cond2.LoverD),(COV_TI1, COV_TI2))
        
    else:
        if COV_WE == None:
            COV_WE = [1 for _ in range(len(z_mesh))]
        else:
            COV_WE = COV_WE * [1 for _ in range(len(z_mesh))]

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
    jgloc           = np.empty(len(z_mesh))
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
        jgloc[0]    = cond.jgloc                                # [m/s]
        jgatm       = cond.jgatm                                # [m/s]

        if preset == 'kim':
            try:
                ai[0]       = cond.area_avg_ai_sheet
                alpha[0]    = cond.area_avg_void_sheet
                Db[0]       = 6 * alpha[0] / ai[0]
            
            except:
                ai[0]       = area_avg(cond,"ai",method=avg_method)                 # [1/m]
                alpha[0]    = area_avg(cond,"alpha",method=avg_method)              # [-]
                Db[0]       = void_area_avg(cond,"Dsm1",method=avg_method) / 1000   # [m]

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
        jgloc[0]    = io["jgloc"][-1]
        jgatm       = io["jgatm"]

    ############################################################################################################################
    #                                                                                                                          #
    #                                                   PRESSURE DROP [Pa/m]                                                   #
    #                                                                                                                          #
    ############################################################################################################################

    # Calculate height change for gravitational loss
    if geometry == 'elbow':
        delta_h = (z_mesh[-1] - z_mesh[0]) * 2 / np.pi          # The height of an elbow is going to be its radius

    elif geometry == 'ubend':
        delta_h = 0

    else:
        delta_h = (z_mesh[-1] - z_mesh[0])                      # Dissipation region is going to be the same as standard VU

    # Calculate initial pressure and pressure gradient
    p = jgatm * p_atm / jgloc[0]                                # Back-calculate local corrected absolute pressure

    if preset == 'kim':
        try:
            p = cond.pz                                         # Override
            dpdz = cond.dpdz

        except:
            dpdz = calc_dpdz(
                cond, 
                method = dpdz_method, 
                chisholm = LM_C, 
                m = m, 
                n = n, 
                k_m = k_m, 
                L = (query - LoverD) * Dh
                ) + ((rho_f * grav * delta_h) / (z_mesh[-1] - z_mesh[0]))

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
    
    jgloc = jgatm * p_atm / pz                                  # Talley used jgP1 and pressure at P1 instead of jgatm and p_atm, but same idea

    ############################################################################################################################
    #                                                                                                                          #
    #                                                           IATE                                                           #
    #                                                                                                                          #
    ############################################################################################################################
    
    # Interpolation-based void fraction estimation methods
    if void_method == 'interp':
        alpha = np.interp(z_mesh / Dh,
            (cond.LoverD, cond2.LoverD),
            (area_avg(cond,'alpha',method=avg_method), area_avg(cond2,'alpha',method=avg_method))
            )

    elif void_method == 'vgz_interp':
        vgz = np.interp(z_mesh / Dh,
            (cond.LoverD, cond2.LoverD),
            (void_area_avg(cond,'ug1',method=avg_method), void_area_avg(cond2,'ug1',method=avg_method))
            )
    
    # Calculate ai(z) to evaluate the steady state one-dim one-group model
    for i, z in enumerate(z_mesh):
        if (i+1) >= len(z_mesh):
            break
        
        ########################################################################################################################
        # Estimate Void Fraction for the current step calculation
        if void_method == 'driftflux':
            pass

        elif void_method == 'continuity':
            pass

        elif void_method == 'interp':
            # alpha[i] already calculated
            pass

        elif void_method == 'vgz_interp':
            # vgz[i] already calculated

            alpha[i] = jgloc[i] / vgz[i]
            Db[i] = 6 * alpha[i] / ai[i]

        elif void_method == 'vgz_talley':
            vgz[i] = 1.05 * (jf + jgloc[i]) - 1.23              # Talley 2012, Eq. 3-31
            alpha[i] = jgloc[i] / vgz[i]
            Db[i] = 6 * alpha[i] / ai[i]

        elif void_method == 'dpdz':
            alpha[i] = calc_void_dpdz(
                cond, jf, Dh, z_mesh, dpdz, LM_C, m, n, k_m, delta_h, grav, mu_f, mu_g, rho_f,
                rho_g = rho_gz[i], L_x = query - LoverD, method = dpdz_method
                )
    
            Db[i] = 6 * alpha[i] / ai[i]

        vgz[i] = jgloc[i] / alpha[i]                            # Estimate void weighted velocity
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
            SEXP[i] = -2 / 3 / pz[i] * ai[i] * vgz[i] * (-dpdz) # Original DOE_MATLAB_IAC
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
        if void_method == 'driftflux':
            j = jgloc[i] + jf
            
            # Drift Velocity
            # Applicable for void fractions less than 20%; for void fractions greater than 30%, use Kataoka and Ishii 1987 for drift-velocity
            vgj = (2**0.5) * (sigma * abs(grav) * (rho_f - rho_gz[i]) / (rho_f**2))**0.25 * (1 - alpha[i])**(1.75)
            
            if C0 == None:
                C0 = C_inf - (C_inf - 1) * np.sqrt(rho_gz[i]/rho_f)     # Round tube drift flux distribution parameter
            
            alpha[i+1] = jgloc[i] / (C0 * j + vgj)

        elif void_method == 'continuity':
            # Original continuity method
            # alpha[i+1] = alpha[i] - alpha[i] / pz[i] * -dpdz * z_step

            # Yadav
            # No discernible difference between this method and original continuity method
            if i <= 2:
                # Specific form of continuity not involving velocity gradients to avoid starting issue
                alpha[i+1] = alpha[i] - alpha[i] * (rho_gz[i+1] - rho_gz[i]) / rho_gz[i]

            else:
                alpha[i+1] = alpha[i] - (alpha[i] / (rho_gz[i] * vgz[i])) * ((rho_gz[i] * vgz[i]) - (rho_gz[i-1] * vgz[i-1]))

        elif void_method == 'interp':
            pass
        
        elif void_method == 'vgz_interp':
            # vgz[i+1] already calculated

            alpha[i+1] = jgloc[i+1] / vgz[i+1]

        elif void_method == 'vgz_talley':
            vgz[i+1] = 1.05 * (jf + jgloc[i+1]) - 1.23               # Talley 2012, Eq. 3-31
            alpha[i+1] = jgloc[i+1] / vgz[i+1]

        elif void_method == 'dpdz':
            alpha[i+1] = calc_void_dpdz(
                cond, jf, Dh, z_mesh, dpdz, LM_C, m, n, k_m, delta_h, grav, mu_f, mu_g, rho_f,
                rho_g = rho_gz[i+1], L_x = query - LoverD, method = dpdz_method
                )

        # Estimate Sauter-mean diameter for the next step calculation
        Db[i+1] = 6 * alpha[i+1] / ai[i+1]
        
    io = {
        "z_mesh"        : z_mesh,
        "ai"            : ai,
        "alpha"         : alpha,
        "pz"            : pz,

        "Db"            : Db,
        "jf"            : jf,
        "jgloc"         : jgloc,
        "jgatm"         : jgatm,
        
        "COV_WE"        : COV_WE,
        "COV_RC"        : COV_RC,
        "COV_TI"        : COV_TI,

        "aiti"          : aiti,
        "airc"          : airc,
        "aiexp"         : aiexp,
        "aiwe"          : aiwe,
        "aivg"          : aivg,
    }

    return io
