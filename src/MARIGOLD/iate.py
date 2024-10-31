from .config import *

def iate_1d_1g(
        # Basic inputs
        cond, query, z_step = 0.01, io = None, geometry = None, cond2 = None,
        
        # IATE Coefficients
        C_WE = None, C_RC = None, C_TI = None, alpha_max = 0.75, C = 3, We_cr = 6, acrit_flag = 0, acrit = 0.13,

        # Method arguments
        preset = None, avg_method = None, cov_method = 'fixed', rf = False, cd_method = 'doe', dpdz_method = 'LM', void_method = 'driftflux',

        # Covariance calculation
        COV_RC = None, COV_TI = None,

        # Pressure drop calculation
        LM_C = 25, k_m = 0.10, m = 0.316, n = 0.25,

        # Void fraction calculation
        C0 = None, C_inf = 1.20,

        ):
    """ Calculate the area-averaged interfacial area concentration at query location based on the 1D 1G IATE
    
    Inputs:
     - cond:            Condition object, part of MARIGOLD framework
     - query:           L/D endpoint
     - z_step:          Axial mesh cell size [-]
     - io:              Output package of iate_1d_1g(), can be used as input for subsequent runs
     - geometry:        Geometry type, defaults to None, can be set to 'elbow', 'ubend'
     - cond2:           Second condition object, for possible interpolation
     - C_WE:            Wake entrainment coefficient
     - C_RC:            Random collision coefficient
     - C_TI:            Turbulent impact coefficient
     - alpha_max:       Maximum void fraction based on HCP bubble distribution, used for random collision calculation
     - C:               C
     - We_cr:           Weber number criterion, used for turbulent impact calculation
     - acrit_flag:      Enable/disable shutting off turbulence-based mechanisms beyond a critical void fraction
     - acrit:           Critical void fraction for shutting off turbulence-based mechanisms
     - preset:          Author preset, fixes coefficients and method arguments to match old MATLAB runs
     - avg_method:      Area-averaging method, can be set to None (for Python Simpson's rule), 'legacy' (for Excel Simpson's Rule)
     - cd_method:       Drag coefficient prediction method, 'err_iter', 'fixed_iter', or 'doe'
     - dpdz_method:     Pressure drop prediction method, 'LM' or 'Kim'
     - void_method:     Void fraction prediction method, 'driftflux' or 'continuity'
     - LM_C:            Lockhart-Martinelli Chisholm parameter
     - k_m:             Minor loss coefficient
     - m:               Friction factor constant
     - n:               Friction factor constant
     - C0:              Drift flux distribution parameter overriding value
     - C_inf:           Drift flux distribution parameter limiting value. Will be used to calculate C0, if none specified

    Notes:
     - Notice some grav terms are made absolute; needs downward flow fixes
     - IATE coefficients set as optional inputs, with default values set depending on geometry
     - vgz calculation in elbow and dissipation length regions still need to be implemented
     - Need a way to compute void fraction across restrictions, void fraction prediction falters
     - Modify MG for Yadav data extraction
     - Revise vgj calculation
    """

    # MARIGOLD retrieval and setup
    theta           = cond.theta                                # Pipe inclination angle [degrees]
    Dh              = cond.Dh                                   # Hydraulic diameter [m]
    rho_f           = cond.rho_f                                # Liquid phase density [kg/m**3]
    rho_g           = cond.rho_g                                # Gas phase density [kg/m**3]
    mu_f            = cond.mu_f                                 # Dynamic viscosity of water [Pa-s]
    mu_g            = cond.mu_g                                 # Dynamic viscosity of air [Pa-s]
    sigma           = cond.sigma                                # Surface tension of air/water [N/m]
    p_atm           = 101325                                    # Ambient pressure [Pa]
    grav            = 9.81*np.sin((theta)*np.pi/180)            # Gravity constant (added by Drew to account for pipe inclination)
    R_spec          = 287.058                                   # Specific gas constant for dry air [J/kg-K]
    T               = 293.15                                    # Ambient absolute temperature [K], for calculating air density as a function of pressure along channel

    # IATE presets (WIP)
    if preset == 'kim':
        theta           = 90
        Dh              = 4*0.20*0.01/2/(0.20+0.01)
        rho_f           = 998
        rho_g           = 1.226
        mu_f            = 0.001
        sigma           = 0.07278
        p_atm           = 101330
        grav            = 9.8

        cd_method       = 'fixed_iter'
        C0              = 1.12
    
    elif preset == 'talley':
        theta           = 0
        Dh              = 0.0381
        rho_f           = 998
        rho_g           = 1.23
        mu_f            = 0.001
        sigma           = 0.07278
        p_atm           = 101353        # Implied by subtracting thesis gauge pressures from MATLAB absolute pressures

        cond.mu_g       = 1.73E-5

        rf              = True
        avg_method      = 'legacy_old'
        cd_method       = 'doe'
        dpdz_method     = 'LM'

        LM_C = 25
    
    elif preset == 'yadav':
        void_method     = 'continuity'

    elif preset == 'worosz':
        cd_method       = 'err_iter'
    
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

    ############################################################################################################################
    #                                                                                                                          #
    #                                                      COEFFICIENTS                                                        #
    #                                                                                                                          #
    ############################################################################################################################

    # For future reference, when adding other geometry types, want to maintain consistency in if/then logic.
    #  1. Check angles with straight pipes
    #  2. Check angles with restrictions
    #  3. Check restrictions
    #       a. 'elbow'
    #       b. 'ubend'
    #       c. 'dissipation'
    #       d. Other
    #  4. Else, default to vertical-upward
    if theta == 0 and geometry == None:         # Horizontal, no elbow (Talley, 2012)
        if C_WE == None:
            C_WE    = 0.000
        if C_RC == None:
            C_RC    = 0.003
        if C_TI == None:
            C_TI    = 0.014
        
        We_cr = 5

    elif geometry == 'elbow':                   # Elbow (Yadav, 2013)
        if C_WE == None:
            C_WE    = 0.000
        if C_RC == None:
            C_RC    = 0.008
        if C_TI == None:
            C_TI    = 0.085

    elif geometry == 'vd':                      # Vertical-downward (Ishii, Paranjape, Kim, and Sun, 2004)
        if C_WE == None:
            C_WE    = 0.002
        if C_RC == None:
            C_RC    = 0.004
        if C_TI == None:
            C_TI    = 0.034

    else:                                       # Default to vertical-upward, no elbow (Ishii, Kim, and Uhle, 2002)
        if C_WE == None:
            C_WE    = 0.002
        if C_RC == None:
            C_RC    = 0.004
        if C_TI == None:
            C_TI    = 0.085
    
    if cond2 == None:
        cov_method = 'fixed'

    if cov_method == 'interp':
        # Use data at initial condition, void reconstruction downstream
        if io == None:
            rf1 = False
            rf2 = rf
        else:
            rf1 = rf
            rf2 = rf
        
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
            # Need to make these into arrays
            # I had logic that did this before... why did I get rid of that
            # Something along the lines of if has property "__len__"?
            COV_RC = 1

        if COV_TI == None:
            COV_TI = 1

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
            ai[0]       = cond.area_avg("ai",method=avg_method)                     # [1/m]
            alpha[0]    = cond.area_avg("alpha",method=avg_method)                  # [-]
            Db[0]       = cond.void_area_avg("Dsm1",method=avg_method) / 1000       # [m]

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
        dpdz = cond.calc_dpdz(
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
                               (cond.void_area_avg('ug1',method=avg_method), cond2.void_area_avg('ug1',method=avg_method))
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
        
        elif geometry == 'elbow':
            SWE[i] = 0
            
            # Yadav calculates vgz differently if solving in elbow region and dissipation length region
            '''
            # Probably want to define these coefficients along with the COEFFICIENT block, if implemented?
            slip = 0.91
            Const1 = 0.05

            beta_diss = 0.18 - 7.6E-7 * Rem     # Dissipation coefficient
            Sratio = 0.1;                       # Sratio = S/S0 (strength of elbow effect)

            # Elbow region
            vgzP4 = slip * jf / (1 - alpha[zstep_v])
            L_D = 0.3156
            deltaz = z[i] - z[zstep_v]
            slope = (vgzP4 - vgz[zstep_v]) / L_D
            vgz[i] = vgz[zstep_v] + slope * (deltaz)
            
            # Dissipation region
            Sratio = np.exp(-beta_diss * (z[i] - z[zstep_v+24]) / Dh)
            vgz[i] = vgzP4 * (1 + Const1 * np.log(Sratio))
            '''

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
                               (cond.void_area_avg('ug1',method=avg_method), cond2.void_area_avg('ug1',method=avg_method))
                               )
            alpha[i+1] = jgloc / vgz[i+1]

        elif void_method == 'pressure_kim':
            f_f, f_g = cond.calc_fric(m = m, n = n)
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
            f_f, f_g = cond.calc_fric(m = m, n = n)
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
