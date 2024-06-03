from .config import *

def iate_1d_1g(
        # Basic inputs
        query, z_step = 0.01,

        # Initial conditions
        cond = None, parcel = None,
        
        # IATE Coefficients
        C_WE = None, C_RC = None, C_TI = None, alpha_max = 0.75, C = 3, We_cr = 6, acrit_flag = 0, acrit = 0.13, C_inf = 1.20,

        # Method arguments
        dpdz_method = 'LM', cd_method = 'doe', void_method = 'driftflux',

        # Pressure drop calculation arguments
        LM_C = 40, k_m = 0.40, LoverD_restriction = None,

        # Temporary arguments
        elbow = False
        ):
    """ Calculate the area-averaged interfacial area concentration at query location based on the 1D 1G IATE

    Version History:
     - v1: Pressure cheating, jgref substitute for jgatm
     - v2: MG update, pressure retrieval, jgatm retrieval
     - v3: Yadav methods, incorporation of COV terms, support for elbows, VU, VD, horizontal
    
    Inputs:
     - cond:            Condition object, part of MARIGOLD framework
     - query:           L/D endpoint
     - z_step:          Axial mesh cell size [-]

     - C_WE:            Wake entrainment coefficient
     - C_RC:            Random collision coefficient
     - C_TI:            Turbulent impact coefficient
     - alpha_max:       Maximum void fraction based on HCP bubble distribution, used for random collision calculation
     - C:               C
     - We_cr:           Weber number criterion, used for turbulent impact calculation
     - acrit_flag:      Enable/disable shutting off turbulence-based mechanisms beyond a critical void fraction
     - acrit:           Critical void fraction for shutting off turbulence-based mechanisms
     - C_inf:           Drift flux distribution parameter limiting value

     - dpdz_method:     Pressure drop prediction method, 'LM' or 'Kim'
     - cd_method:       Drag coefficient prediction method, 'iter' or 'doe'
     - void_method:     Void fraction prediction method, 'driftflux' or 'continuity'

    Notes:
     - Isn't this script nice and clean? Follow good coding practices, kids. - David
     - IATE coefficients set as optional inputs, with default values set depending on geometry
     - COV terms being implemented in Condition.py, not incorporated into this function yet
     - vgz calculation in elbow and dissipation length regions still need to be implemented
     - There must be something wrong with drift flux void fraction calculation
     - Modify MG for Yadav data extraction
    """

    # MARIGOLD retrieval
    if cond != None and parcel == None:
        theta       = cond.theta                                # Pipe inclination angle
        Dh          = cond.Dh                                   # Hydraulic diameter
        LoverD      = cond.LoverD                               # Condition L/D

        rho_f       = cond.rho_f                                # Liquid phase density [kg/m**3]
        rho_g       = cond.rho_g                                # Gas phase density [kg/m**3]
        mu_f        = cond.mu_f                                 # Viscosity of water [Pa-s]
        sigma       = cond.sigma                                # Surface tension of air/water [N/m]
    else:
        theta       = parcel["theta"]
        Dh          = parcel["Dh"]
        LoverD      = parcel["z_mesh"][-1] / Dh

        rho_f       = parcel["rho_f"]
        rho_g       = parcel["rho_g"]
        mu_f        = parcel["mu_f"]
        sigma       = parcel["sigma"]
    
    p_atm           = 101325                                    # Ambient pressure [Pa]
    grav            = 9.81*np.sin((theta)*np.pi/180)            # Gravity constant (added by Drew to account for pipe inclination)

    R_spec          = 287.058                                   # Specific gas constant for dry air [J/kg-K]
    T               = 293.15                                    # Ambient absolute temperature [K], for calculating air density as a function of pressure along channel
    
    # Mesh generation
    if query < LoverD+z_step:
        raise ValueError('Please choose a query L/D downstream of the boundary condition.')
    
    z_mesh = np.arange(LoverD, query + z_step, z_step)          # Axial mesh [-]
    z_mesh = z_mesh * Dh                                        # Axial mesh [m], units necessary for dp calculation
    z_step = z_step * Dh

    # Yadav
    '''
    RoverD          = 9                                         # Make into optional input
    R_curv          = RoverD * Dh                               # Radius of curvature of elbow/bend
    L_curv          = np.pi/2 * R_curv                          # This assumes 90\degree elbow    

    # So, Yadav defines under Expcond.m the length of the vertical section -- I think this is how he toggles between VU, H, VD, elbow
    # Instead, we probably would want to have an input switch indicating which mode of IATE we're executing
    '''

    ############################################################################################################################
    #                                                                                                                          #
    #                                                      COEFFICIENTS                                                        #
    #                                                                                                                          #
    ############################################################################################################################

    if theta == 90 and elbow == False:                          # Vertical, no elbow
        if C_WE == None:
            C_WE    = 0.002
        if C_RC == None:
            C_RC    = 0.004        
        if C_TI == None:
            C_TI    = 0.085

        COV_RC      = 1
        COV_TI      = 1

    elif theta == 0 and elbow == False:                         # Horizontal, no elbow, look at Ran's work?
        if C_WE == None:
            C_WE    = 0.000
        if C_RC == None:
            C_RC    = 0.003
        if C_TI == None:
            C_TI    = 0.014
        
        COV_RC      = 1
        COV_TI      = 1

    elif elbow == True:                                         # Elbow, look at Shoxu's work?
        if C_WE == None:
            C_WE    = 0.000
        if C_RC == None:
            C_RC    = 0.012
        if C_TI == None:
            C_TI    = 0.085
        
        COV_RC      = 1
        COV_TI      = 1

    else:
        if C_WE == None:
            C_WE    = 0.002
        if C_RC == None:
            C_RC    = 0.004        
        if C_TI == None:
            C_TI    = 0.085
        
        COV_RC      = 1
        COV_TI      = 1
    
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

    if cond != None and parcel == None:
        aiwe[0]     = 0
        airc[0]     = 0
        aiti[0]     = 0
        aiexp[0]    = 0
        aivg[0]     = 0

        ai[0]       = cond.area_avg("ai")                       # [1/m]
        alpha[0]    = cond.area_avg("alpha")                    # [-]
        Db[0]       = cond.void_area_avg("Dsm1") / 1000         # [m]
        jf          = cond.jf                                   # [m/s]
        jgloc       = cond.jgloc                                # [m/s]
        jgatm       = cond.jgatm                                # [m/s]

        # Pressure drop [Pa/m]
        if LoverD_restriction == None:
            L = None
        else:
            L = LoverD_restriction*Dh

        dpdz        = cond.calc_dpdz(method=dpdz_method, LM_C=LM_C, k_m=k_m, L=L)
    else:
        aiwe[0]     = parcel["aiwe"][-1]
        airc[0]     = parcel["airc"][-1]
        aiti[0]     = parcel["aiti"][-1]
        aiexp[0]    = parcel["aiexp"][-1]
        aivg[0]     = parcel["aivg"][-1]

        ai[0]       = parcel["ai"][-1]
        alpha[0]    = parcel["alpha"][-1]
        Db[0]       = parcel["Db"][-1]
        jf          = parcel["jf"]
        jgloc       = parcel["jgloc"]
        jgatm       = parcel["jgatm"]

        dpdz        = parcel["dpdz"]
    
    # Local Pressure along the test section
    p = (jgatm * p_atm / jgloc) - p_atm                         # Back-calculate local corrected gauge pressure
    pz = (p + p_atm) * (1 - (z_mesh - z_mesh[0]) * abs(dpdz / (p + p_atm)))
    
	# Local gas density along the test section
    rho_gz = rho_g * pz / p_atm                                 # Talley
    # rho_gz = pz / R_spec / T                                  # Worosz, Ideal Gas Law
    
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
        vgz[i] = jgloc / alpha[i]                               # Estimate void weighted velocity
        vfz = jf / (1 - alpha[i])

        ########################################################################################################################
        # Estimate bubble relative velocity <ur> (See Talley, 2012, 4.2.2.6)
        ur = 0.2                                                    # Set to constant value of 0.23 m/s by Schilling (2007)

        if cd_method == 'iter':
            err = 0.1
            while abs(err) > 0.000001:
                ReD = rho_f * ur * Db[i] * (1 - alpha[i]) / mu_f    # Bubble Reynolds number
                CDe = 24 * (1 + 0.1 * ReD**0.75) / ReD              # Drag coefficient

                # Relative velocity (Ishii and Chawla, 1979), implemented by Worosz accounting for density difference
                ure1 = (4 * grav * Db[i] * (rho_f - rho_gz[i]) * (1 - alpha[i]) / 3 / CDe / rho_f)**0.5
                err = (ure1 - ur) / ur

                ur = ure1
            ReD = rho_f * ur * Db[i] * (1 - alpha[i]) / mu_f        # Update bubble Reynolds number
            CDwe = 24 * (1 + 0.1 * ReD**0.75) / ReD                 # Update drag coefficient

        elif cd_method == 'doe':
            # Original DOE_MATLAB_IAC
            ReD = rho_f * ur * Db[i] * (1 - alpha[i]) / mu_f
            CDwe = 24 * (1 + 0.1 * ReD**0.75) / ReD
            ur = (4 * grav * Db[i] / 3 / CDwe)**0.5                 # Interestingly, Yadav keeps 9.8 instead of changing grav for angle

        ########################################################################################################################
        # Estimate Energy Dissipation Rate and Turbulent Velocity (See Talley, 2012, 4.2.2.3)
        #   > One-group models written using turbulent fluctuation velocity, while models implemented in TRACE are written using
        #     dissipation rate

        mu_m = mu_f / (1 - alpha[i])                            # Mixture viscosity, given by Ishii and Chawla (Eq. 4-10 in Kim, 1999)
        rho_m = (1 - alpha[i]) * rho_f + alpha[i] * rho_gz[i]   # Mixture density

        vm = (rho_f * (1 - alpha[i]) * vfz + rho_gz[i] * alpha[i] * vgz[i]) \
            / rho_m                                             # Mixture velocity
        
        # Rem = rho_f * vm * Dh / mu_m                          # Mixture Reynolds number
        Rem = rho_m * vm * Dh / mu_m                            # Yadav
        fTW = 0.316 * (mu_m / mu_f)**0.25 / Rem**0.25           # Two-phase friction factor
        e = fTW * (vm**3) / 2 / Dh                              # Energy dissipation rate (Wu et al., 1998; Kim, 1999)
        ut = 1.4 * e**(1/3) * Db[i]**(1/3)                      # Turbulent velocity (Batchelor, 1951; Rotta, 1972)

        ########################################################################################################################
        # Estimate sources & sinks in the Interfacial Area Transport Eqn. (Part 1)

        # Sink due to Wake Entrainment
        if theta == 90 and elbow == False:            
            SWE[i] = C_WE * CDwe**(1/3) * ur * ai[i]**2 / 3 / np.pi

        elif theta == 0 and elbow == False:
            SWE[i] = 0
        
        elif elbow == True:
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
        
        # Sink due to Random Collisions	
        RC1 = ut * ai[i]**2 / alpha_max**(1/3) / (alpha_max**(1/3) - alpha[i]**(1/3))
        RC2 = 1 - np.exp(-C * alpha_max**(1/3) * alpha[i]**(1/3) / (alpha_max**(1/3) - alpha[i]**(1/3)))
        SRC[i] = COV_RC * C_RC * RC1 * RC2 / 3 / np.pi
        
        # Source due to Turbulent Impact
        TI1 = ut * ai[i]**2 / alpha[i]
        We = rho_f * ut**2 * Db[i] / sigma                      # Weber number criterion

        if We > We_cr:
            TI2 = (1 - We_cr / We)**0.5 * np.exp(-We_cr / We)
        else:
            TI2 = 0
        STI[i] = COV_TI * C_TI * TI1 * TI2 / 18

        ########################################################################################################################
        # Estimate sources & sinks in the Interfacial Area Transport Eqn. (Part 2)

        # Source due to Bubble Expansion
        if i <= 2:      # Previously 3, but in MATLAB (1 indexing vs. 0 indexing)
            # Forward difference for first node
            SEXP[i] = -2 / 3 / rho_gz[i] * ai[i] * vgz[i] * (rho_gz[i+1] - rho_gz[i]) / z_step
        else:
            # Backwards difference for remaining nodes
            SEXP[i] = -2 / 3 / rho_gz[i] * ai[i] * vgz[i] * (rho_gz[i] - rho_gz[i-1]) / z_step

        # SEXP[i] = -2 / 3 / pz[i] * ai[i] * vgz[i] * (-dpdz)   # Original DOE_MATLAB_IAC
        
        # Source/sink due to Bubble Acceleration (advection in Yadav's script)
        if i <= 2:
            dvg = 0
            dvgdz = 0

            # dvgdz = dvg / z_step                              # Gonna leave this here for a good chuckle
        else:
            dvg = vgz[i] - vgz[i-1]
            dvgdz = dvg / z_step

        SVG[i] = ai[i] * dvgdz

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

        # Estimate Void Fraction for the next step calculation

        if void_method == 'driftflux':      # Drift Flux Model

            j = jgloc + jf
            
            # Drift Velocity
            # Applicable for void fractions less than 20%; for void fractions greater than 30%, use Kataoka and Ishii 1987 for drift-velocity
            vgj = (2**0.5) * (sigma * grav * (rho_f - rho_gz[i]) / (rho_f**2))**0.25 * (1 - alpha[i])**(1.75)
            
            C0 = C_inf - (C_inf - 1) * np.sqrt(rho_gz[i]/rho_f)     # Round tube drift flux distribution parameter
            # alpha[i+1] = (jgloc) / (C0 * j + vgj)

            alpha[i+1] = (jgloc) / (C0 * j + vgj)

        elif void_method == 'continuity':   # Continuity

            # alpha[i+1] = alpha[i] - alpha[i] / pz[i] * -dpdz * z_step

            # Yadav
            # No discernible difference between this method and original continuity method
            if i <= 2:
                # Specific form of continuity not involving velocity gradients to avoid starting issue
                alpha[i+1] = alpha[i] - alpha[i] * (rho_gz[i+1] - rho_gz[i]) / rho_gz[i]

            else:
                alpha[i+1] = alpha[i] - (alpha[i] / (rho_gz[i] * vgz[i])) * ((rho_gz[i] * vgz[i]) - (rho_gz[i-1] * vgz[i-1]))

        # Estimate Sauter mean diameter for the next step calculation
        Db[i+1] = 6 * alpha[i+1] / ai[i+1]
        
    parcel = {
        "theta"         : theta,
        "Dh"            : Dh,
        "LoverD"        : LoverD,
        "rho_f"         : rho_f,
        "rho_g"         : rho_g,
        "mu_f"          : mu_f,
        "sigma"         : sigma,
        "aiti"          : aiti,
        "airc"          : airc,
        "aiexp"         : aiexp,
        "aiwe"          : aiwe,
        "aivg"          : aivg,
        "ai"            : ai,
        "alpha"         : alpha,
        "Db"            : Db,
        "jf"            : jf,
        "jgloc"         : jgloc,
        "jgatm"         : jgatm,
        "dpdz"          : dpdz,

        "z_mesh"        : z_mesh,
        "pz"            : pz
    }

    return parcel