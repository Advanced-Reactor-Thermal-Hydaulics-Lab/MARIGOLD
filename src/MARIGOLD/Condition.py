from .config import *
from .plot_utils import *
from scipy import interpolate
from scipy.optimize import minimize


class Condition:
    """ Class to handle the local probe data

    Data is stored in the Condition.data property. It's actually 3 layers of dictionary
    data [angle] gives a dictionary with the various r/R
    data [angle][r/R] gives a dictionary with the MIDAS output
    The MIDAS output is itself a dictionary, with the keys listed in the "tab_keys" array
    So data[angle][r/R]['alpha'] should give you the void fraction at r/R for phi = angle
    This structure is initialized with zeros for the MIDAS output at the pipe center and wall

    Can also get the data at a local point from calling the condition, syntax

    cond(phi, r, 'param'). Phi is in radians, the arguments can be constants or numpy arrays.
    Also has an option for interpolation, 'interp_method' 

"""

    debugFID = None
    def __init__(self, jgref:float, jgloc:float, jf:float, theta:int, port:str, database:str, fluids = 'air-water', g = 9.81) -> None:
        """ Initialize Condition object

        Inputs:
        - jgref, reference superficial gas velocity
        - jgloc, local superficial gas velocity
        - jf, superficial liquid velocity
        - theta, angle of inclination of flow direction (0° is horizontal, 90° is vertical upwards)
        - port, string to denote the port
        - database, string to associate what database this data is from

        Optional inputs:
        - fluids, what fluid pair to use as the gas and liquid
        - g, gravitational acceleration. In case you're on Mars. MARIGOLD multiplies by sin(θ)

        Implemented Fluids:
        - air-water, uses properties at atmospheric conditions
        
        """
        self.jgref = jgref
        self.jf = jf
        self.jgloc =jgloc
        self.theta = theta
        self.port = port
        self.database = database

         # Check if port is one of the specified vertical downward ports after U-bend, set theta to -90 degree if true (Quan 10/25)
        # if self.port in {"P5", "P6", "P7", "P8", "P9", "P10"}:
        #    self.theta = int(-90)
        # else:
        #     self.theta = int(90)

        self.name = f"jf={self.jf}_jgloc={self.jgloc:0.2f}_theta={self.theta}_port={self.port}_{self.database}"

        # Data is stored in this phi array. 3 layers of dictionary
        # phi [angle] gives a dictionary with the various r/R
        # phi [angle][r/R] gives a dictionary with the MIDAS output
        # So phi[angle][r/R]['alpha'] should give you the void fraction at r/R for phi = angle
        # This structure is initialized with zeros for the MIDAS output at the pipe center and wall
        self._angles = np.arange(0, 361, 22.5) # HARDCODED 22.5 degree increments
        #self.phi = deepcopy(dict( zip(angles, deepcopy([ {0.0: dict( zip(tab_keys, [0]*len(tab_keys)) ), 1.0: dict(zip(tab_keys, [0]*len(tab_keys)) ) } ]) * len(angles)) ))
        self.data = {}

        self.mirrored = False
        self.FR = 0 # Flow regime variable. 0 is undefined, 1 is bubbly, etc.

        self.j = self.jgloc + self.jf

        if 'D' in self.port and 'P' not in self.port:
            self.LoverD = int(self.port.strip('D'))
        else: # Assume it's PITA
            if self.port == 'P1':
                self.LoverD = 30
            elif self.port == 'P2':
                self.LoverD = 66
            elif self.port == 'P3':
                self.LoverD = 110
            elif self.port == 'P4':
                self.LoverD = 130.04
            elif self.port == 'P5':
                self.LoverD = 144.17
            elif self.port == 'P6':
                self.LoverD = 147.57
            elif self.port == 'P7':
                self.LoverD = 151.84
            elif self.port == 'P8':
                self.LoverD = 165.57
            elif self.port == 'P9':
                self.LoverD = 194.07          
            elif self.port == 'P10':
                self.LoverD = 230.07
            else:
                self.LoverD = -1
                print(f"Warning: Could not determine port L/D for {self}")

        self.vwvg = -1
        self.void_cov = -1

        self.area_avg_void_sheet = -1

        if database == 'Ryan' or database == 'ryan' or database == 'adix':
            self.Dh = 0.0254 # m, for Ryan
            self.marker_type = 'o'
            self.marker_color = 'r'
            self.line_style = 'solid'

        elif database == 'Kong':
            self.Dh = 0.1016
            self.marker_type = 's'
            self.marker_color = 'b'
            self.line_style = 'dotted'

        elif database == 'Talley':
            self.Dh = 0.0381
            self.marker_type = '^'
            self.marker_color = 'g'
            self.line_style = 'dashed'

        elif database == 'Yadav':
            self.Dh = 0.0508
            self.marker_type = 'D'
            self.marker_color = 'purple'
            self.line_style = 'dashdot'

        else:
            print(f"Warning: Could not determine Dh for {self.name}")
            self.Dh = np.nan
            self.marker_type = '$?$'
            self.marker_color = 'yellow'

        # Empty dictionaries, filled when max or area avg is called
        self.area_avgs = {}
        self.circ_seg_area_avgs = {}
        self.maxs = {}
        self.mins = {}
        self._grads_calced = []

        if fluids == 'air-water' or fluids == 'water-air':
            self.rho_f = 998       # kg/m^3
            self.rho_g = 1.204     # kg/m^3
            self.mu_f = 0.001002   # Pa s
            self.mu_g = 0.01803e-3 # Pa s

            self.sigma = 0.0728    # N/m

            self.Ref = self.rho_f * self.jf * self.Dh / self.mu_f
            self.Reg = self.rho_g * self.jgloc * self.Dh / self.mu_g
        
        else:
            raise NotImplementedError(f"{fluids} not available, try 'air-water'")
        
        self.g = g
        self.gz = g * np.sin(self.theta * np.pi/180)
        return

    def __eq__(self, __o: object) -> bool:
        if isinstance(__o, Condition):

            return ((self.jf == __o.jf) and (self.jgref == __o.jgref) and (self.theta == __o.theta) and (self.port == __o.port) and self.database == __o.database)
        
        return False

    def __hash__(self) -> int:
        return hash(repr(self))

    def __repr__(self) -> str:
        return self.name

    def __call__(self, phi_in:np.ndarray, r_in:np.ndarray, param:str, interp_method='None') -> np.ndarray:
        """Returns the value of param at (phi, r). Phi is in radians, r nondimensional
         
        Interp method:
         - 'None', will try to fetch raw data at this location
         - 'linear', linear interpolation
         - 'spline', spline interpolation
         - 'linear_xy', cartesian interpolation, phi -> x, r -> y
           
        """
        if type(phi_in) != np.ndarray:
            if debug: warnings.warn("Converting phi_in to np.ndarray")
            phi_in = np.asarray(phi_in)

        if type(r_in) != np.ndarray:
            if debug: warnings.warn("Converting r_in to np.ndarray")
            r_in = np.asarray(r_in)

        
        if interp_method == 'None':
            try:
                param_values = np.zeros((r_in.size, phi_in.size)) # To-do: check if r_in and phi_in actually exist
                for i, r_val in enumerate(r_in):
                    for j, phi_val in enumerate(phi_in):
                        try:
                            param_values[i,j] = self.data[round(float(phi_val) * 180 / np.pi, 2)][r_val][param]
                        except KeyError as e:
                            if abs(abs(r_val) - 1) < 0.0001:
                                param_values[i,j] = 0
                            else:
                                print(e)
                                raise 
            except:
                # Probably input a single phi instead of an array
                
                param_values = self.data[round(float(phi_in) * 180 / np.pi, 2)][float(r_in)][param]
            return param_values
        
        elif interp_method == 'spline':
            try:
                return self.spline_interp(phi_in, r_in)
            except:
                self.mirror(uniform_rmesh=True)
                self.fit_spline(param)
                return self.spline_interp[param](phi_in, r_in)
        
        elif interp_method == 'linear':
            try:
                return self.linear_interp[param](phi_in, r_in)
            except:
                self.fit_linear_interp(param)
                return self.linear_interp[param](phi_in, r_in)
            
        elif interp_method == 'linear_xy':
            x = phi_in
            y = r_in
            try:
                return self.linear_xy_interp[param](x, y)
            except:
                self.fit_linear_xy_interp(param)
                return self.linear_xy_interp[param](x, y)
        
        else:
            raise NameError(f"{interp_method} not recognized. Accepted arguments are 'None', 'spline', 'linear' or 'linear_xy'")

    def pretty_print(self, print_to_file= False, FID=debugFID, mirror=False) -> None:
        """Prints out all the information in a Condition in a structured way 

        Specifically, everything in the Condition.phi dictionary, which has angles
        and r/Rs.

        Can either print to a file (specified by FID) or to stdout. Option to mirror the 
        data, if that hasn't already been done

        """

        print(f"jf = {self.jf}\tjg = {self.jgref}\ttheta = {self.theta}\t{self.port}\t{self.database}", file=FID)
        
        if mirror:
            self.mirror()
        
        if print_to_file:
            for angle, r_dict in self.data.items():
                print(angle, file=FID)
                for r, midas_output in r_dict.items():
                    print(f"\t{r}", file=FID)
                    print("\t\t", midas_output, file=FID)
        else:
            for angle, r_dict in self.data.items():
                print(angle)
                for r, midas_output in r_dict.items():
                    print(f"\t{r}")
                    print("\t\t", midas_output)
        return

    def mirror(self, method = 'sym90', sym90 = False, axisym = False, uniform_rmesh = False, uniform_rmesh_fill = 'interp', force_remirror=False) -> None:
        """ Mirrors data, so we have data for every angle

        :param method: method to use for mirroring. Options are 'sym90', 'axisym', and 'avg_axisym'. Defaults to 'sym90'
        :type method: str, optional
        :param sym90: overwrites method, for backwards compatibility, defaults to False
        :type sym90: bool, optional
        :param axisym: overwrites method, for backwards compatibility, defaults to False
        :type axisym: bool, optional
        :param uniform_rmesh: ensure every angle has data for every r/R point. Will linearly interpolate when data on either side is available. This only  considers the +r/R mesh, defaults to False
        :type uniform_rmesh: bool, optional
        :param force_remirror: Usually can only mirror once. This overrides that restriction, defaults to False
        :type force_remirror: bool, optional

        Saves:
         - original_mesh
        """

        # Only ever call this function once
        if self.mirrored and not force_remirror:
            return
        
        # Have a record of actual measurement locations

        # First step is to find the phi angles that have data        
        angles_with_data = set()
        self.original_mesh = []

        for angle, rdict in self.data.items():
            for rstar, midas_dict in rdict.items():
                if any(midas_dict.values()):
                    if 'num_spherical' in midas_dict.keys():
                        if (abs(midas_dict['num_spherical'] - 1894) < 0.01) and (abs(midas_dict['num_cap'] - 11) < 0.01):
                            # this is dummy data, ignore
                            continue
                        if (abs(midas_dict['num_spherical'] - 302) < 0.01) and (abs(midas_dict['ai_distorted'] - 1.25) < 0.01):
                            # this is dummy data, ignore
                            continue
                        if (abs(midas_dict['num_spherical'] - 157) < 0.01) and (abs(midas_dict['ai_distorted'] - 0.88) < 0.01):
                            # this is dummy data, ignore
                            continue

                    angles_with_data.add(angle)
                    self.original_mesh.append( (angle, rstar) )

        angles_with_data = list(angles_with_data)
        # print(angles_with_data)
        if debug: print('Angles with data: ', angles_with_data, file=debugFID)
        num_angles_measured = len(angles_with_data) # the actual number of angles measured (r/R can be negative)

        # Next, take any negative values for the existing data and copy them to the 
        # +180 phi angle (the complementary angle)
        # Also delete the negative data from the existing phi, so all angles only have postive entries
        angles_to_add = []
        all_rs = set([0.0])
        for angle in angles_with_data:
            if angle <= 180:
                comp_angle = angle + 180

            if (comp_angle not in angles_with_data) and (comp_angle <= 360) and (comp_angle not in angles_to_add):

                data = deepcopy(self.data[angle])
                rs = list(deepcopy(self.data[angle]).keys())

                for r in rs:
                    
                    if r > 0:
                        all_rs.add(r)
                        data.pop(r)
                
                for r in rs:
                    if r < 0:
                        data[-r] = self.data[angle].pop(r)
                        data.pop(r)
                    
                    elif r == 0:
                        pass
                
                self.data[angle].update({1.0: deepcopy(zero_data)}) 

                # There should always be data at r/R 0 so we can plot contours
                try: 
                    dummy = self.data[angle][0.0]
                except:
                    self.data[angle].update({0.0: deepcopy(zero_data)})  # just in case


                self.data.update({comp_angle: {}})
                self.data[comp_angle].update( {1.0: deepcopy(zero_data)} )
                self.data[comp_angle].update( data )
                

                angles_to_add.append(comp_angle)

        angles_with_data += angles_to_add
        if debug: print('Angles with data after comp_angle: ', angles_with_data, file=debugFID)
        
        if (360 not in (angles_with_data)) and (0 in angles_with_data):
            ref_angle = 0
            data = deepcopy(self.data[ref_angle])
            self.data.update({360: {}})
            self.data[360].update( data )
            angles_with_data.append(360)

        # Now comes the actual mirroring step. Need data for every angle, incremements of 22.5° (self._angles)
        if sym90:
            method == 'sym90'
        elif axisym:
            method == 'axisym'

        if method == 'axisym': 
            # axisymmetric
            # Find the reference angle with the most data
            ref_angle = angles_with_data[0] # initial guess
            for angle in angles_with_data:
                if len(self.data[ref_angle].keys()) < len(self.data[angle].keys()):
                    ref_angle = angle
            
            for angle in self._angles:
                data = deepcopy( self.data[ref_angle] )
                
                data.update({1.0: deepcopy(zero_data)}) 
                self.data.update({angle: {}})
                self.data[angle].update( data )

        elif method.lower() == 'avg_axisym':
            # print('avg_axisym')
            # Find a reference angle with a complement, replace all the data with the average of the two
            ref_angle = angles_with_data[0] # initial guess
            for angle in angles_with_data:
                if (angle+180 in angles_with_data) and len(self.data[ref_angle].keys()) < len(self.data[angle].keys()):
                    ref_angle = angle

            temp_data = deepcopy(self.data[ref_angle])
            comp_data = deepcopy(self.data[ref_angle+180])
            avg_data = deepcopy(self.data[ref_angle])
            # print('comp data', comp_data)
            for rstar, midas_dict in temp_data.items():
                avg_data.update({rstar: {} })
                for param, val in midas_dict.items():
                    try:
                        avg_data[rstar][param] = (val + comp_data[rstar][param]) / 2
                        # print('used avg data')
                    except:
                        avg_data[rstar][param] = val

            for angle in self._angles:
                data = deepcopy(avg_data)
                self.data.update({angle: {}})
                self.data[angle].update( data )

        elif method.lower() == 'sym90': 
            # symmetric across the 90 degree line
            for angle in self._angles:
                if angle not in angles_with_data:

                    data = {0.0: deepcopy(zero_data), 1.0: deepcopy(zero_data)} # Fine if this gets overwritten, just need to make sure there's some data at 0 for plotting

                    if angle <= 90:
                        # Quadrant I, should usually have data here, but if we don't, try to copy data from Q2
                        ref_angle = 180 - angle
                        try:
                            data = deepcopy(self.data[ref_angle])
                        except KeyError:
                            if debug: print(f"No data found for {angle} when mirroring {self.name}, defaulting to 0s")
                            data = {0.0: deepcopy(zero_data)}

                    elif angle > 90 and angle <= 180:
                       # Quadrant II, mirror from Quadrant I
                        ref_angle = 180 - angle
                        data = deepcopy(self.data[ref_angle])

                    elif angle > 180 and angle <= 270:
                        # Quadrant III, should be covered by the negative of Quadrant I
                        # But if we're here there's no data here. So make sure it's 0
                        ref_angle = 540 - angle
                        try:
                            data = deepcopy(self.data[ref_angle])
                        except KeyError:
                            if debug: print(f"No data found for {angle} when mirroring {self.name}, defaulting to 0s")
                            data = {0.0: deepcopy(zero_data)}

                    elif angle > 270 and angle < 360:
                        # Quadrant IV, mirror from Quadrant III
                        ref_angle = 540 - angle
                        data = deepcopy(self.data[ref_angle])

                    elif angle == 360:
                        ref_angle = 0
                        data = deepcopy(self.data[ref_angle])

                    data.update({1.0: deepcopy(zero_data)}) # paranoia
                    
                    # Check if data exists at zero, and if not, just put some zero data in
                    try: 
                        dummy = data[0.0]
                    except:
                        data.update({0.0: deepcopy(zero_data)}) # just in case
                    
                    if angle > 360: continue # Just in case
                    self.data.update({angle: {}})
                    self.data[angle].update( data )

            # Check if data exists at zero, and if not, just put some zero data in
                try: 
                    dummy = self.data[angle][0.0]
                except:
                    self.data[angle].update({0.0: deepcopy(zero_data)}) # just in case


        else:
            # No symmetry being assumed. But we still want data at every angle. If it doesn't exist, must be 0
            warnings.warn("Filling in all angles without data as 0")
            for angle in self._angles:
                if angle not in angles_with_data:
                    data = {0.0: deepcopy(zero_data)}
                    data.update({1.0: deepcopy(zero_data)})
                    if angle > 360: continue
                    self.data.update({angle: {}})
                    self.data[angle].update( data )

        if uniform_rmesh:
            for angle in self.data.keys():
                if angle not in self._angles:
                    self._angles.append(angle)
            for param in tab_keys:
                self.fit_linear_interp(param)
            #print(all_rs)
            # Make sure all the angles have data for all the rpoints
            self.all_rs = list(all_rs)
            self.all_rs.sort()
            #print(self.all_rs)
            for angle in self._angles:
                for r in all_rs:
                    if r not in self.data[angle].keys(): # This will break if there's data for 0.85 in some cases but not others
                        self.data[angle].update({r: deepcopy(zero_data)})
            
            # so go through and interpolate the points where we have data on either side
            # for angle in self._angles:
            #     for i in range(len(self.all_rs) - 2):
            #         #print(angle, self.all_rs[i+1])
            #         if (self.data[angle][self.all_rs[i+2]]['alpha'] != 0) and (self.data[angle][self.all_rs[i]]['alpha'] != 0) and (self.data[angle][self.all_rs[i+1]]['alpha'] == 0):
            #             print(f"Warning: interpolating data for {angle}°, {self.all_rs[i+1]} to maintain uniform r/R mesh")
            #             for param in tab_keys:
            #                 try:

            #                     x = self.all_rs[i+1]
            #                     x1 = self.all_rs[i+2]
            #                     y1 = self.data[angle][x1][param]
            #                     x2 = self.all_rs[i]
            #                     y2 = self.data[angle][x2][param]

            #                     interp = y1 + (y2-y1)/(x2-x1) * (x - x1)
            #                     self.data[angle][self.all_rs[i+1]][param] = interp
                            
            #                 except KeyError:
            #                     #print(f"{param} not found for {angle}°, {self.all_rs[i+1]}")
            #                     continue
                        
            #             self.data[angle][self.all_rs[i+1]]['roverR'] = f"interpolated, {angle}, {i+1}"
            #         else:
            #             pass
            #             #print(f"I'm not interpolating for {angle}°, {self.all_rs[i+1]}")


            for phi in self._angles:
                for r in all_rs:
                    for param in tab_keys:
                        try:
                            if self.data[phi][r][param] == 0 and r < 1.0:
                                if (uniform_rmesh_fill) == 0:
                                    self.data[phi][r][param] = 0
                                elif uniform_rmesh_fill == 'interp':
                                    self.data[phi][r][param] = self(phi*np.pi/180, r, param, interp_method = 'linear')
                        except Exception as e:
                            if debug: print("Error interpolating in unifrom rmesh, ", phi, r, e)
                            pass

            for param in tab_keys:
                self.fit_linear_interp(param) # Because it was created without the filled-in data before

        self.mirrored = True

        # clean up nones
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                for param, value in midas_dict.items():
                    if value is None:
                        print(f"Warning: 'None' value for {param} at {rstar}. Setting to 0")
                        midas_dict.update( {param : 0})

        return
    
    def add_mesh_points(self, r_points:list, suppress=False):
        """Method for adding additional r/R points

        Data linearly interpolated based on surrounding data (specifically using __call__ at the (angle, r) location in question)
        
        Inputs:
         - r_points, points to add
        
        """

        for angle in self.data.keys():
            for r in r_points:
                temp_midas_data = []
                
                r_keys = list(self.data[angle].keys())
                r_keys.sort()

                params = self.data[angle][r_keys[-2]]
                
                for param in params:
                    try:
                        temp_midas_data.append( float(self(angle*np.pi/180, r, param, interp_method = 'linear' ) ))
                    except KeyError as e:
                        if suppress == False:
                            print(e)
                        
                        temp_midas_data.append( 0 )
                
                self.data[angle].update( {r: dict(zip(tab_keys, temp_midas_data))} )

        return 1
    
    def approx_vf(self, n=7, overwrite_vf = False) -> None:
        """Method for approximating vf with power-law relation. 

        .. math:: v_{f, approx} = \\frac{(n+1)(2*n+1)}{ (2*n^{2})} * (j_{f} / (1- \\langle \\alpha \\rangle)) * (1 - abs(rstar))**(1/n)

        Stores:
         - "vf_approx" in midas_dict
         - "vf", if it does not alread exist in midas_dict

        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                vf_approx = (n+1)*(2*n+1) / (2*n*n) * (self.jf / (1-self.area_avg('alpha'))) * (1 - abs(rstar))**(1/n)
                if 'vf' in midas_dict.keys():
                    if debug: print(f"approx_vf: data found for {angle}\t{rstar}", file=debugFID)
                    if overwrite_vf:
                        midas_dict.update({'vf': vf_approx})
                else:
                    midas_dict.update({'vf': vf_approx})
                
                midas_dict.update({'vf_approx': vf_approx})

        return
    
    def approx_vg(self, method = 'vr', n=7) -> None:
        """Method for approximating vg with power-law relation. I don't think this makes sense

        .. math:: v_{g, approx} = \\frac{(n+1)(2*n+1)}{ (2*n^{2})} * (j_{g} / \\langle \\alpha \\rangle) * (1 - abs(rstar))**(1/n)

        Stores:
         - "vg_approx" in midas_dict
         - "ug1", if it does not alread exist in midas_dict

        If you have a lot of group II bubbles this functions' no good

        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                if method == 'power-law':
                    vg_approx = (n+1)*(2*n+1) / (2*n*n) * (self.jgloc / (self.area_avg('alpha'))) * (1 - abs(rstar))**(1/n)
                elif method == 'vrmodel':
                    if 'vr_model' not in midas_dict.keys():
                        self.calc_vr_model()
                    vg_approx = midas_dict['vf_approx'] + midas_dict['vr_model']
                
                if 'ug1' in midas_dict.keys():
                    if debug: print(f"approx_vg: data found for {angle}\t{rstar}", file=debugFID)
                    vg_approx = midas_dict['ug1']
                else:
                    midas_dict.update({'ug1': vg_approx})
                
                midas_dict.update({'vg_approx': vg_approx})

        return
    
    def approx_vf_Kong(self, n=7) -> None:
        """Method for approximating vf from Kong. TODO 

        Not currently implemented

        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                vf_approx = (n+1)*(2*n+1) / (2*n*n) * (self.jf / (1-self.area_avg('alpha'))) * (1 - abs(rstar))**(1/n)
                midas_dict.update({'vf': vf_approx})

        return
    
    def calc_vf_lee(self, K=1):
        """ Calculate vf, jf, vr based on Lee et al. (2002) equation

        Inputs:
         - None

        Stores:
         - "vf_lee" in midas_dict
         - "jf_lee" in midas_dict
         - "vr_lee" in midas_dict

        Returns:
         - Area-average vf_lee

        Really a model proposed by Bosio and Malnes

        .. math:: v_{f} = \\frac{ 1 }{ \\sqrt{1 - \\alpha^{2} / 2} } * \\sqrt{ \\frac{ 2 \\Delta p }{K \\rho_{f}} } 
        
        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    dp = midas_dict['delta_p'] * 6894.757
                except:
                    raise NotImplementedError("Δp needed for cacluclation of vf_lee")
                
                vf_lee = 1 / np.sqrt(1 - midas_dict['alpha']**2/2) * np.sqrt( 2 * dp / (K * self.rho_f))
                
                midas_dict.update({'vf_lee': vf_lee})
                midas_dict.update({'jf_lee': (1-midas_dict['alpha'])* vf_lee})
                
                vg = midas_dict['ug1']
                if vg == 0:
                    vr_lee = 0
                else:
                    vr_lee = vg - vf_lee
                midas_dict.update({'vr_lee':  vr_lee})

        return self.area_avg('vf_lee')
    
    def calc_vf_naive(self):
        """ Calculate vf, jf, vr based on single-phase equation

        Inputs:
         - None

        Stores:
         - "vf_naive" in midas_dict
         - "jf_naive" in midas_dict
         - "vr_naive" in midas_dict

        Returns:
         - Area-average vf_naive

        .. math::  v_{f} = \\sqrt{ \\frac{ 2 \\Delta p }{ \\rho_{f}} } 
        
        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    vf_naive = midas_dict['vf_naive']
                except KeyError:
                    try:
                        dp = midas_dict['delta_p'] * 6894.757
                    except:
                        raise NotImplementedError("Δp needed for cacluclation of vf_naive")
                    
                    vf_naive = np.sqrt( 2*dp / self.rho_f)
                
                midas_dict.update({'vf_naive': vf_naive})
                midas_dict.update({'jf_naive': (1-midas_dict['alpha'])* vf_naive})
                
                vg = midas_dict['ug1']
                if vg == 0:
                    vr_naive = 0
                else:
                    vr_naive = vg - vf_naive
                midas_dict.update({'vr_naive':  vr_naive})

        return self.area_avg('vf_naive')
    
    def calc_vr(self, quiet = False) -> None:
        """Method for calculating relative velocity. 
        
        Inputs:
         - warn_approx, flag for warning if :math`v_r` is not found and function is approximating

        Note that if vg = 0, then this method says vr = 0. This will happen when no data is present,
        such as in the bottom of the pipe in horizontal, when this is not necessarily true

        Stores:
         - 'vr' in midas_dict

        Returns:
         - Area-average vr

        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    vf = midas_dict['vf']
                except:
                    if not quiet:
                        print("Warning: Approximating vf in calculating vr, since no data found")
                        warn_approx = False
                    self.approx_vf()
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
                    

        return self.area_avg('vr')
    
    def calc_vr2(self, warn_approx = True) -> None:
        """Method for calculating relative velocity based on ug2. 
        
        Inputs:
         - warn_approx, flag for warning if :math`v_r` is not found and function is approximating

        Note that if vg2 = 0, then this method says vr2 = 0. This will happen when no data is present,
        such as in the bottom of the pipe in horizontal, when this is not necessarily true

        Stores:
         - 'vr2' in midas_dict

        Returns:
         - Area-average vr2

        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    vf = midas_dict['vf']
                except:
                    if warn_approx:
                        print("Warning: Approximating vf in calculating vr, since no data found")
                        warn_approx = False
                    self.approx_vf()
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
                    

        return self.area_avg('vr2')

    def calc_vgj(self, warn_approx = True) -> None:
        """Method for calculating Vgj

        Inputs:
         - warn_approx, flag for warning if :math`V_{gj}` is not found and function is approximating

        Stores:
         - "vgj" in midas_dict
         - "j" in midas_dict
         - "alpha_j" in midas_dict

        Returns:
         - Void-weighted area-averaged Vgj

        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    dummy = midas_dict['vf']
                except:
                    if warn_approx:
                        print("Warning: Approximating vf in calculating local j, since no data found")
                        warn_approx = False
                    self.approx_vf()
                
                j_local = midas_dict['alpha'] * midas_dict['ug1'] + (1 - midas_dict['alpha']) * midas_dict['vf']
                vgj = midas_dict['ug1'] - j_local
                alpha_j = midas_dict['alpha'] * j_local
                midas_dict.update({'vgj': vgj})
                midas_dict.update({'j': j_local})
                midas_dict.update({'alpha_j': alpha_j})

        return self.void_area_avg('vgj')

    def calc_grad(self, param: str, recalc = False) -> None:
        """Calculates gradient of param based on the data in self. 

        Inputs:
         - "param", string for local parameter to calculate the gradient of
        
        Stores:
         - grad\_'param'_r
         - grad\_'param'_phi
         - grad\_'param'_phinor
         - grad\_'param'_x
         - grad\_'param'_y
         - grad\_'param'_total

        Returns:
         - Area average grad\_'param'_total

        """
        
        if not self.mirrored: self.mirror()

        if param in self._grads_calced and not recalc:
            return
        else:
            self._grads_calced.append(param)

        grad_param_name = 'grad_' + param

        phis = self._angles.tolist()
        phis.sort()
        maxj = len(phis)

        for phi_angle, r_dict in self.data.items():
            rs = list(r_dict.keys())

            rs.sort()
            #print(rs)

            for i in range(1, len(rs) ):
                if rs[i] < 0: 
                    if debug: print(f'Warning: somehow we still have negative data.\n{self}\n{phi_angle=}\t{rs=}')
                grad_r_param = (r_dict[rs[i]][param] - r_dict[rs[i-1]][param]) / (rs[i] - rs[i-1])
                
                
                j = phis.index(phi_angle)
                if j < maxj-1:
                    hi_phi = phis[(j+1)]
                else:
                    hi_phi = phis[1] # use 22.5
                
                lo_phi = phis[(j-1) ]

                if j == 0:
                    lo_phi = phis[(j-2) ] # use 347.5, not 360
                    #print(hi_phi, lo_phi)

                # print(phis[j], hi_phi, lo_phi)
                
                try:
                    hi = self.data[hi_phi][rs[i]][param]
                except KeyError as e:
                    if debug: print(f"Key error found when indexing {e} for hi. Likely a case of the data being zero for the adjacent point, setting to 0...", file=debugFID)
                    hi = 0

                try:
                    lo = self.data[lo_phi][rs[i]][param]
                except KeyError as e:
                    if debug: print(f"Key error found when indexing {e} for lo. Likely a case of the data being zero for the adjacent point, setting to 0...", file=debugFID)
                    lo = 0

                if rs[i]> 0:
                    grad_phi_param = 1./rs[i] * (hi - lo) / ((hi_phi - lo_phi) * np.pi/180)
                else:
                    grad_phi_param = 0 # I guess? Shouldn't actually come up

                r_dict[rs[i]].update( {grad_param_name+'_r': grad_r_param } )
                r_dict[rs[i]].update( {grad_param_name+'_phi': grad_phi_param } )
                r_dict[rs[i]].update( {grad_param_name+'_phinor': grad_phi_param*rs[i] } )

                if rs[i] != 0:
                    r_dict[rs[i]].update( {grad_param_name+'_y': grad_r_param * np.sin(phi_angle*180/np.pi) + np.cos(phi_angle*180/np.pi)/(rs[i])*grad_phi_param } )
                    r_dict[rs[i]].update( {grad_param_name+'_x': grad_r_param * np.cos(phi_angle*180/np.pi) - np.sin(phi_angle*180/np.pi)/(rs[i])*grad_phi_param } )

                r_dict[rs[i]].update( {grad_param_name+'_total': grad_r_param+grad_phi_param } )
                
                
                if i == 1: # first point, also calculate derivative at 0, using grad_r_param
                
                    # x and y only a function of r, dψ/dφ = 0
                    if phi_angle == 0:
                        # r_dict[0.0].update( {grad_param_name+'_x': grad_r_param * np.cos(phi_angle) } )
                        # Copy to all other angles
                        for temp_angle in self._angles:
                            self.data[temp_angle][0.0].update({grad_param_name+'_x': grad_r_param * np.cos(phi_angle*180/np.pi) } )
                    
                    elif phi_angle == 90:
                        # r_dict[0.0].update( {grad_param_name+'_y': grad_r_param * np.sin(phi_angle) } )
                        # Copy to all other angles
                        for temp_angle in self._angles:
                            self.data[temp_angle][0.0].update({grad_param_name+'_y': grad_r_param * np.sin(phi_angle*180/np.pi) } )


                    if phi_angle <= 180:
                        comp_angle = phi_angle + 180
                    else:
                        comp_angle = phi_angle - 180
                    
                    # Want to average dψ/dr between current phi angle and complementary angle

                    if grad_param_name+'_r' in r_dict[0.0].keys():
                        value = deepcopy(r_dict[0.0][grad_param_name+'_r']) + 0.5 * grad_r_param
                        r_dict[0.0].update( {grad_param_name+'_r': value} )
                    else:
                        r_dict[0.0].update( {grad_param_name+'_r': 0.5 * grad_r_param} )

                    
                    if grad_param_name+'_r' in self.data[comp_angle][0.0].keys():
                        value = deepcopy(self.data[comp_angle][0.0][grad_param_name+'_r']) + 0.5 * grad_r_param
                        self.data[comp_angle][0.0].update( {grad_param_name+'_r': value} )
                    else:
                        self.data[comp_angle][0.0].update( {grad_param_name+'_r': 0.5 * grad_r_param} )
                        
                    
                    r_dict[0.0].update( {grad_param_name+'_phi': 0 } )     
                    self.data[comp_angle][0.0].update({grad_param_name+'_phi': 0 })

                    r_dict[0.0].update( {grad_param_name+'_phinor': 0 } )     
                    self.data[comp_angle][0.0].update({grad_param_name+'_phinor': 0 })

                    r_dict[0.0].update( {grad_param_name+'_total': deepcopy(r_dict[0.0][grad_param_name+'_r'])} )
                    self.data[comp_angle][0.0].update({grad_param_name+'_total': deepcopy(r_dict[0.0][grad_param_name+'_r']) })
                    

        return self.area_avg(grad_param_name+'_total')
    
    def fit_spline(self, param: str) -> None:
        """Fits a RectBivariateSpline for the given param. 

           Can access later with self.spline_interp[param]. Must specify the 'param' to fit
        
        """
        
        if not self.check_param(param, strict=False):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")


        try: dummy = self.spline_interp
        except:
            self.spline_interp = {}
        
        if param in self.spline_interp.keys():
            if debug: print(f"{param} already has a fit spline")
            return

        self.mirror(uniform_rmesh=True)
        rs = []
        phis = []
        vals = []

        phi_unique = set()
        r_unique = set()

        for angle, r_dict in self.data.items():
            phi_unique.add(angle)
            for rstar, midas_dict in r_dict.items(): # TODO check for the 0.85 issue
                r_unique.add(rstar)

        r_unique = np.asarray(list(r_unique))
        phi_unique = np.asarray(list(phi_unique))
        r_unique.sort()
        phi_unique.sort()
        
        for angle in phi_unique:
            for rstar in r_unique:
                #print(angle, rstar, self.phi[angle][rstar][param])
                vals.append(self.data[angle][rstar][param])

        #print(r_unique, phi_unique)

        Vals = np.asarray(vals).reshape((phi_unique.size, r_unique.size))
        #print(Vals)

        #print(phi_knots, r_knots)
        spline_interpolant = interpolate.RectBivariateSpline(phi_unique* np.pi/180, r_unique, Vals, kx=3)
        self.spline_interp.update({param: spline_interpolant})
        return
    
    def fit_linear_interp(self, param: str, suppress=True, refit = False) -> None:
        """Makes a LinearNDInterpolator for the given param. 
        
            Access with self.linear_interp[param], phi in radians

            If the requested param does not exist for a given (r, phi), this fills in 0
             
        """

        if not self.check_param(param, strict=False):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")


        try: dummy = self.linear_interp
        except:
            self.linear_interp = {}
        
        if param in self.linear_interp.keys() and not refit:
            if debug: print(f"{param} already has a linear interpolator")
            return
        
        self.mirror()

        rs = []
        phis = []
        vals = []

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                rs.append(rstar)
                phis.append(angle * np.pi/180)
                try:
                    vals.append(float(midas_dict[param]))
                except KeyError:
                    
                    if not suppress:
                        print(self.name, angle, rstar, "has no ", param)
                        raise(KeyError)
                    
                    vals.append( 0 )        # Temporary, figure this out later

        linear_interpolant = interpolate.LinearNDInterpolator(list(zip(phis, rs)), vals, fill_value=0)
        self.linear_interp.update({param: linear_interpolant})

        return
    
    def fit_linear_xy_interp(self, param: str, suppress=True) -> None:
        """ Makes a LinearNDInterpolator for the given param in x y coords
        
        Can access  with self.linear_xy_interp[param]
                          
        """

        if not self.check_param(param):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")

        try: dummy = self.linear_xy_interp
        except:
            self.linear_xy_interp = {}
        
        if param in self.linear_xy_interp.keys():
            if debug: print(f"{param} already has a linear interpolator in xy space")
            return
        
        self.mirror()

        xs = []
        ys = []
        vals = []

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                xs.append(rstar * np.cos(angle * np.pi / 180))
                ys.append(rstar * np.sin(angle * np.pi / 180))
                try:
                    vals.append(midas_dict[param])
                except KeyError:
                    if suppress == False:
                        print(self.name, angle, rstar, "has no ", param)

        linear_interpolant = interpolate.LinearNDInterpolator(list(zip(xs, ys)), vals, fill_value=0)
        self.linear_xy_interp.update({param: linear_interpolant})

        return
    
    def sum(self, param: str) -> float:
        """ Adds all values of param
        
        Inputs:
         - param, string of local parameter to sum

        Returns:
         - sum of the value of parameter at every point
                          
        """
        total = 0
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                total += midas_dict[param]

        return total
    
    def max(self, param: str, recalc=False) -> float:
        """ Return maximum value of param in the Condition
        
        Inputs:
         - param, string of local parameter to find maximum of
         - recalc, if already calculated, will pull from dictionary instead of calculating unless this is true

        Returns:
         - maximum value of param
                          
        """
        
        if (param in self.maxs.keys()) and (not recalc):
            return self.maxs[param] # why waste time 
        max = 0
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] > max:
                    max = midas_dict[param]
                    location = rstar
        self.maxs.update({param:max})
        return (max)

    def max_loc(self, param: str)-> tuple:
        """ Return location of maximum value of param
        
        Inputs:
         - param, string of local parameter to find maximum of

        Returns:
         - r/R location of the maximum value of param
                          
        """
        max = 0
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] > max:
                    max = midas_dict[param]
                    location = (rstar, angle)

        return (location)
    
    def min(self, param: str, recalc = True, nonzero=False)-> float:
        """ Return minimum value of param in the Condition
        
        Inputs:
         - param, string of local parameter to find minimum of
         - recalc, if already calculated, will pull from dictionary instead of calculating unless this is true
         - nonzero, if true, will exclude points with the param = 0 when finding the minimum

        Returns:
         - minimum value of param
                          
        """

        if (param in self.mins.keys()) and (not recalc):
            return self.mins[param] # why waste time 
        min = 10**7
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                
                if nonzero and midas_dict[param] < 1e-7:
                    continue
                
                if midas_dict[param] < min:
                    min = midas_dict[param]
                    location = rstar
        self.mins.update({param:min})
        return (min)

    def min_loc(self, param: str)-> float:
        """ Return location of minimum value of param
        
        Inputs:
         - param, string of local parameter to find minimum of

        Returns:
         - r/R location of the maximum value of param
                          
        """
        min = 10**7
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] < min:
                    min = midas_dict[param]
                    location = (rstar, angle)

        return (location)
    
    def max_line_loc(self, param: str, angle:float) -> float:
        """ Return location of maximum value of param at a given angle
        
        Inputs:
         - param, string of local parameter to find maximum of
         - angle, angle to search along

        Returns:
         - r/R location of the maximum value of param when φ=angle
                          
        """
        max = 0
        for rstar, midas_dict in self.data[angle].items():
            if midas_dict[param] > max:
                max = midas_dict[param]
                location = rstar

        return (location)
    
    def max_line(self, param: str, angle:float) -> float:
        """ Return maximum value of param at a given angle
        
        Inputs:
         - param, string of local parameter to find maximum of
         - angle, angle to search along

        Returns:
         - maximum value of param when φ=angle
                          
        """
        max = 0
        for rstar, midas_dict in self.data[angle].items():
            if midas_dict[param] > max:
                max = midas_dict[param]
                location = rstar

        return (max)
    
    def find_hstar_pos(self, method='max_dsm', void_criteria = 0.05) -> float:
        """ Returns the vertical distance from the top of the pipe to the bubble layer interface

        Inputs:
         - method
         - void_criteria, used in some methods
        
        Methods for determining bubble layer interface
         - max_dsm
         - min_grad_y
         - max_grad_y
         - max_mag_grad_y
         - zero_void
         - percent_void, search down the phi = 90 line, find largest value of rstar where alpha < min_void
         - Ryan_Ref, uses 1.3 - 1.57e-5 * Ref, proposed by Ryan (2022)

        Returns:
         - hstar
        
        """

        if method == 'max_dsm':
            r_max, phi_max = self.max_loc('Dsm1')
            h_star= 1 - r_max * np.sin(phi_max*np.pi/180)
            return h_star
        
        elif method == 'min_grad_y':
            self.calc_grad('alpha')
            r_min, phi_min = self.min_loc('grad_alpha_y')
            h_star= 1 - r_min * np.sin(phi_min*np.pi/180)
            return h_star

        elif method == 'max_grad_y':
            self.calc_grad('alpha')
            r_min, phi_min = self.max_loc('grad_alpha_y')
            h_star= 1 - r_min * np.sin(phi_min*np.pi/180)
            return h_star
        
        elif method == 'max_mag_grad_y':
            self.calc_grad('alpha')
            for rstar, midas_dict in self.data[90].items():
                if midas_dict['alpha'] < void_criteria:
                    pass
                    #TODO finish this
            h_star= 1 - roverRend
            return h_star

        elif method == 'zero_void':
            # search down the phi = 90 line, find largest value of rstar where alpha < min_void
            roverRend = -1
            max_loc = self.max_line_loc('alpha', 90)
            max_void = self.max_line('alpha', 90)
            for rstar, midas_dict in self.data[90].items():
                if midas_dict['alpha'] < void_criteria:
                    if rstar > roverRend and rstar < max_loc:
                        roverRend = rstar
            
            if roverRend == -1: # Didn't find it on the positive side
                for rstar, midas_dict in self.data[270].items():
                    if midas_dict['alpha'] < void_criteria:
                        if -rstar > roverRend:
                            roverRend = -rstar

            self.roverRend = roverRend
            h_star = 1 - self.roverRend
            return h_star
        
        elif method == 'percent_void':
            # search down the phi = 90 line, find largest value of rstar where alpha < min_void
            roverRend = -1
            max_loc = self.max_line_loc('alpha', 90)
            max_void = self.max_line('alpha', 90)
            for rstar, midas_dict in self.data[90].items():
                if midas_dict['alpha'] < max_void * void_criteria:
                    if rstar > roverRend and rstar < max_loc:
                        roverRend = rstar
            
            if roverRend == -1:
                for rstar, midas_dict in self.data[180].items():
                    if midas_dict['alpha'] < max_void * void_criteria:
                        if -rstar > roverRend:
                            roverRend = -rstar

            self.roverRend = roverRend
            h_star = 1 - self.roverRend
            return h_star
        
        elif method == 'Ryan_Ref':
                        
            self.roverRend = 1.3 - 1.57e-5 * self.Ref
            h_star = 1 - self.roverRend
            return h_star

        print('Invalid method for find_h_pos')
        return np.NaN

    def avg(self, param: str, include_zero=True) -> float:
        """Calculates a basic average of a parameter, "param"
        
        Inputs:
         - param, string of local parameter to find average of
         - include_zero, if you want to include when param = 0 in average

        Returns:
         - average value of param

        """
        count = 0
        avg_param = 0
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                
                if include_zero:
                    count += 1
                    avg_param += midas_dict[param]
                else:
                    if abs(midas_dict[param]) > 0:
                        count += 1
                        avg_param += midas_dict[param]

        return avg_param / count
    
    def check_param(self, param:str, strict=True) -> bool:
        """Checks if a parameter is present in condition data. Ignores r/R = 1.0

        :param param: parameter to check
        :type param: str
        :param strict: Check each point. Otherwise, one value every angle is acceptable., defaults to True
        :type strict: bool, optional
        :return: True if param is present, false if not
        :rtype: bool
        """        
        points_to_add = []
        for angle, r_dict in self.data.items():
            param_in_angle = False
            
            for rstar, midas_dict in r_dict.items():
                if param not in midas_dict.keys():
                    if strict:
                        if abs(rstar - 1.0) < 0.001:
                            points_to_add.append( (angle, rstar))
                            pass # TODO add extra stuff if it doesn't exist at 1.0, it really shouldn't be a problem
                        else:
                            # return (False, rstar, angle)
                            return False
                
                else:
                    param_in_angle = True

            if (not param_in_angle) and (not strict):
                # return (False, rstar, angle)
                return False

        return True
    
    def check_param_loc(self, param:str) -> tuple:
        """ Checks if a given parameter is present everywhere in midas_dict
        
        """

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                if param not in midas_dict.keys():
                    return (angle, rstar)

        return ()

    def area_avg(self, param: str, even_opt='first', recalc=True, method=None) -> float:
        """Method for calculating the area-average of a parameter, "param"
        
        Inputs:
         - param, string of local parameter to area-average
         - even_opt, passed to integration method
         - recalc, if area-average previously calculated, it was stored in a dictionary, and can just pull from there

        Returns:
         - Area-averaged parameter, or None if the method failed
        
        Uses Simpson's rule for integration, specifically scipy.integrate.simpsons
        
        """

        # Check that the parameter that the user requested exists
        if not self.check_param(param):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")
        
        if (param in self.area_avgs.keys()) and (not recalc):
            return self.area_avgs[param] # why waste time, if we already calculated this don't do it again
        
        if method == 'legacy':
            I = 0
            param_r = [] # array for parameter integrated wrt r
            angles = []
            
            if not self.mirrored:
                warnings.warn("Mirroring in area-avg")
                self.mirror()

            for angle, r_dict in self.data.items():

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
            self.area_avgs.update({param: I})

        elif method == 'legacy_old':
            # There's definitely a more elegant way to do this, but I just want Talley's data to work for now

            I = 0
            param_r = [] # array for parameter integrated wrt r
            angles = []
            
            if not self.mirrored:
                warnings.warn("Mirroring in area-avg")
                self.mirror()

            for angle, r_dict in self.data.items():
                if angle == 360:
                    continue

                try:
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

                except Exception as e:
                    print(e)

            I = sum(param_r) / 8
            self.area_avgs.update({param: I})

        else:
            # We have to integrate twice, once with resepect to r, again with respect to phi
            # Start with r
            I = 0
            param_r = [] # array for parameter integrated wrt r
            angles = []
            
            if not self.mirrored:
                warnings.warn("Mirroring in area-avg")
                self.mirror()

            for angle, r_dict in self.data.items():

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
            self.area_avgs.update({param: I})

        return I

    def circ_segment_area_avg(self, param:str, hstar:float, ngridr=25, ngridphi=25, int_err = 10**-4) -> float:
        """Method for calculating the area-average of a parameter, "param" in some circular segment defined by :math:`h^*`
        
        Inputs:
         - param, string of local parameter to area-average
         - hstar, distance from the top of the pipe that defines circular segment
         - int_err, acceptable warning for integral error (passed to scipy.integrate.dblquad)

        Returns:
         - Area-averaged parameter, or None if the method failed
        
        Uses scipy.integrate.dblquad. May be computationally expensive
        
        """

        # Check that the parameter that the user requested exists
        # Check that the parameter that the user requested exists
        if not self.check_param(param):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")
        
        self.mirror()

        rs = []
        phis = []
        vals = []
        for angle, r_dict in self.data.items():
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

    def circ_segment_void_area_avg(self, param:str, hstar:float, ngridr=25, ngridphi=25, int_err = 10**-4) -> float:
        """Method for calculating the void-weighted area-average of a parameter, "param" in some circular segment defined by :math:`h^*`
        
        Inputs:
         - param, string of local parameter to area-average
         - hstar, distance from the top of the pipe that defines circular segment
         - int_err, acceptable warning for integral error (passed to scipy.integrate.dblquad)

        Returns:
         - Void-weighted rea-averaged parameter, or None if the method failed
        
        Uses scipy.integrate.dblquad. May be computationally expensive
        
        """

        # Check that the parameter that the user requested exists
        if not self.check_param(param):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")
        
        
        self.mirror()

        rs = []
        phis = []
        vals = []
        denom = []
        for angle, r_dict in self.data.items():
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

    def line_avg(self, param:str, phi_angle:float, even_opt='first') -> float:
        """Method for calculating the average value of param across a diameter defined by φ = angle
        
        Inputs:
         - param, string of local parameter to area-average
         - phi_angle, angle to define diameter to integrate across
         - even_opt, see :any:area_avg

        Returns:
         - Average value of param along φ = angle
        
        Uses scipy.integrate.simpsons
        
        """

        # Check that the parameter that the user requested exists
        self.mirror()

        if phi_angle not in self.data.keys():
            if debug: print(self.data, file=debugFID)
            print(f"Could not area-average {param} for condition {self.name}\nData for {phi_angle} not found after mirroring!")
            return


        if not self.check_param(param):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")

        r_for_int = []
        var_for_int = []

        for rstar, midas_dict in self.data[phi_angle].items():
            if rstar not in r_for_int:
                r_for_int.append(rstar)
                var_for_int.append(midas_dict[param])

        if phi_angle <=180:
            comp_angle = phi_angle+180
            
        else:
            comp_angle = phi_angle - 180
        
        for rstar, midas_dict in self.data[comp_angle].items():
            if rstar not in r_for_int:
                r_for_int.append(-rstar)
                var_for_int.append(midas_dict[param])

        var_for_int = [param for _, param in sorted(zip(r_for_int, var_for_int))]
        r_for_int = sorted(r_for_int)

        I = integrate.simpson(y=var_for_int, x=r_for_int) / 2 # Integrate wrt theta, divide by normalized length

        return I

    def line_avg_dev(self, param:str, phi_angle:float, even_opt='first') -> float:
        """Method for calculating the average value of param across a diameter defined by φ = angle
        
        Inputs:
         - param, string of local parameter to calculate line average deviation
         - phi_angle, angle to define diameter to integrate across
         - even_opt, see :any:area_avg

        Returns:
         - integrand result
        
        Taking param = :math:`\\psi`, mathematically, this function does

        .. math:: \\frac{\\langle \\psi - \\langle \\psi^2 \\rangle  \\rangle }{\\langle \\psi \\rangle^2}
        
        """

        # Check that the parameter that the user requested exists
        self.mirror()

        if phi_angle not in self.data.keys():
            if debug: print(self.data, file=debugFID)
            print(f"Could not area-average {param} for condition {self.name}\nData for {phi_angle} not found after mirroring!")
            return


        if not self.check_param(param):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")


        r_for_int = []
        var_for_int = []

        for rstar, midas_dict in self.data[phi_angle].items():
            if rstar not in r_for_int:
                r_for_int.append(rstar)
                var_for_int.append((midas_dict[param] - self.area_avg(param))**2)

        if phi_angle <=180:
            comp_angle = phi_angle+180
            
        else:
            comp_angle = phi_angle - 180
        
        for rstar, midas_dict in self.data[comp_angle].items():
            if rstar not in r_for_int:
                r_for_int.append(-rstar)
                var_for_int.append((midas_dict[param] - self.area_avg(param))**2)

        var_for_int = [param for _, param in sorted(zip(r_for_int, var_for_int))]
        r_for_int = sorted(r_for_int)

        I = integrate.simpson(y=var_for_int, x=r_for_int) / 2 / self.area_avg(param)**2 # Integrate wrt theta, divide by normalized length

        return I


    def void_area_avg(self, param: str, even_opt='first', method = None) -> float:
        """Method for calculating the void-weighted area-average of a parameter, "param"
        
        Inputs:
         - param, string of local parameter to area-average
         - even_opt, passed to integration method
         - recalc, if area-average previously calculated, it was stored in a dictionary, and can just pull from there

        Returns:
         - void-weighted area-averaged parameter, or None if the method failed
        
        Uses Simpson's rule for integration, specifically scipy.integrate.simpsons
        
        """

        # Check that the parameter that the user requested exists
        if not self.check_param(param):
            raise KeyError(f"Invalid parameter {param} selected. Not present at {self.check_param_loc(param)}")


        if method == 'legacy':
            I = 0
            param_r = [] # array for parameter integrated wrt r
            angles = []
            
            if not self.mirrored:
                warnings.warn("Mirroring in area-avg")
                self.mirror()

            for angle, r_dict in self.data.items():

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

            I = sum(param_r) / 8 / self.area_avg('alpha',method=method)

        elif method == 'legacy_old':
            # We have to integrate twice, once with resepect to r, again with respect to phi
            # Start with r

            I = 0
            param_r = [] # array for parameter integrated wrt r
            angles = []
            
            if not self.mirrored:
                warnings.warn("Mirroring in area-avg")
                self.mirror()

            for angle, r_dict in self.data.items():
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
                    
            I = sum(param_r) / 8 / self.area_avg('alpha',method=method)

        else:
            # We have to integrate twice, once with resepect to r, again with respect to phi
            # Start with r

            I = 0
            param_r = [] # array for parameter integrated wrt r
            angles = []
            
            self.mirror()


            for angle, r_dict in self.data.items():

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

            I = integrate.simpson(y=param_r, x=angles) / np.pi / self.area_avg('alpha') # Integrate wrt theta, divide by normalized area

        return I
    
    def interp_area_avg(self, param:str, interp_type = 'linear', int_error = 10**-6) -> float:
        """Method for calculating the area-average of a parameter, "param", based on an interpolation of the data
        
        Inputs:
         - param, string of local parameter to area-average
         - interp_type, option for what kind of interpolation to use

        Options for interp_type:
         - linear, see :any:`fit_linear_interp`
         - spline, see :any:`fit_spline`

        Returns:
         - area-averaged parameter
        
        Uses scipy.integrate.dblquad, may be computationally expensive
        
        """

        def integrand(phi, r): # phi will be in radians from dblquad
            return self(phi, r, param, interp_method=interp_type) * r
        
        I = integrate.nquad(integrand, ranges = [(0, 1), (0, np.pi * 2)], opts = {'epsabs': int_error, 'points': list(self.data[0].keys()), 'points': list(self.data.keys())})[0] / self.interp_area_avg('alpha') / np.pi
        return I
    
    def interp_void_area_avg(self, param:str, interp_type = 'linear', int_error = 10**-6) -> float:
        """Method for calculating the void-weighted area-average of a parameter, "param", based on an interpolation of the data
        
        Inputs:
         - param, string of local parameter to area-average
         - interp_type, option for what kind of interpolation to use

        Options for interp_type:
         - linear, see :any:`fit_linear_interp`
         - spline, see :any:`fit_spline`

        Returns:
         - void weighted area-averaged parameter
        
        Uses scipy.integrate.dblquad, may be computationally expensive
        
        """

        def integrand(phi, r): # phi will be in radians from dblquad
            return self(phi, r, param, interp_method=interp_type) * self(phi, r, 'alpha', interp_method=interp_type) * r
        
        I = integrate.nquad(integrand, ranges = [(0, 1), (0, np.pi * 2)], opts = {'epsabs': int_error, 'points': list(self.data.keys()), 'points': list(self.data[0].keys())} )[0] / self.interp_area_avg('alpha') / np.pi
        return I

    def spline_void_area_avg(self, param:str, int_error = 10**-6) -> float:
        """Function to void-weighted area-average param based on a spline interpolation

        Inputs:
         - param, string of local parameter to area-average
        
        """

        def integrand(phi, r):
            return self.spline_interp[param](phi * 180/np.pi, r) * self.spline_interp['alpha'](phi * 180/np.pi, r) * r
        
        def integrand_denom(phi, r):
            return self.spline_interp['alpha'](phi * 180/np.pi, r) * r
        
        I = integrate.dblquad(integrand, 0, 1, 0, np.pi * 2, epsabs = int_error)[0] / integrate.dblquad(integrand_denom, 0, 1, 0, np.pi * 2, epsabs = int_error)[0]
        return I
    
    def spline_circ_seg_area_avg(self, param:str, hstar:float, int_err = 10**-4) -> float:
        """Function to area-average over a circular segment using the spline interpolation of param

        Inputs:
         - param, string of local parameter to area-average
         - hstar, distance from the top of the pipe that defines circular segment
         - int_err, acceptable warning for integral error (passed to scipy.integrate.dblquad)

        Returns:
         - integrand result

        Uses scipy.integrate.dblquad. May be computationally expensive
        
        """

        def integrand(r, phi):
            return self.spline_interp[param](phi * 180/np.pi, r) * r

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

    def calc_void_cov(self):
        """Calculates the void covariance

        Inputs:
         - None

        Stores:
         - "void_cov" in self

        Returns:
         - void covariance

        Mathematically performing the operation
        
        .. math:: \\frac{\\langle \\alpha^2 \\rangle }{\\langle \\alpha \\rangle^2}
        
        """

        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        for angle, r_dict in self.data.items():
            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                    vars_temp.append(rstar * midas_dict['alpha']**2)
                    if debug: print(angle, midas_dict, file=debugFID)
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            if debug: print("Arrays to integrate", rs, vars, file=debugFID)
                
            param_r.append( integrate.simpson(y=vars, x=rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(y=vars, x=rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)
        param_r_int = [var for _, var in sorted(zip(angles, param_r))]
        angles_int = sorted(angles)
        I = integrate.simpson(y=param_r_int, x=angles_int) / np.pi / self.area_avg('alpha')**2 # Integrate wrt theta, divide by normalized area

        self.void_cov = I
        return I

    def calc_sigma_alpha(self):
        """Calculates the second moment of alpha

        Inputs:
         - None

        Stores:
         - "sigma_alpha" in self

        Returns:
         - second moment of :math:`alpha`

        Mathematically performing the operation

        .. math:: \\frac{\\langle (\\alpha - \\langle \\alpha \\rangle)^2 \\rangle }{\\langle \\alpha \\rangle^2}
        
        """

        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        alpha_avg = self.area_avg('alpha')

        for angle, r_dict in self.data.items():
            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                    vars_temp.append(rstar * (midas_dict['alpha'] - alpha_avg)**2)
                    if debug: print(angle, midas_dict, file=debugFID)
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            if debug: print("Arrays to integrate", rs, vars, file=debugFID)
                
            param_r.append( integrate.simpson(y=vars, x=rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(y=vars, x=rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)
        param_r_int = [var for _, var in sorted(zip(angles, param_r))]
        angles_int = sorted(angles)
        I = integrate.simpson(y=param_r_int, x=angles_int) / np.pi / alpha_avg**2 # Integrate wrt theta, divide by normalized area
        if debug: print('Calculated sigma_alpha: ', I)

        self.sigma_alpha = I
        return I
    
    def calc_sigma_alpha1(self):
        """Calculates the second moment of alpha_G1

        Inputs:
         - None

        Stores:
         - "sigma_alpha1" in self

        Returns:
         - second moment of :math:`alpha_G1`

        Mathematically performing the operation

        .. math:: \\frac{\\langle (\\alpha - \\langle \\alpha \\rangle)^2 \\rangle }{\\langle \\alpha \\rangle^2}
        
        """

        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        alpha_avg = self.area_avg('alpha_G1')

        for angle, r_dict in self.data.items():
            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                    vars_temp.append(rstar * (midas_dict['alpha_G1'] - alpha_avg)**2)
                    if debug: print(angle, midas_dict, file=debugFID)
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            if debug: print("Arrays to integrate", rs, vars, file=debugFID)
                
            param_r.append( integrate.simpson(vars, rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(vars, rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)
        param_r_int = [var for _, var in sorted(zip(angles, param_r))]
        angles_int = sorted(angles)
        I = integrate.simpson(param_r_int, angles_int) / np.pi / alpha_avg**2 # Integrate wrt theta, divide by normalized area
        if debug: print('Calculated sigma_alpha1: ', I)

        self.sigma_alpha1 = I
        return I
    
    def calc_sigma_ug(self):
        """Calculates the second moment of ug1

        Inputs:
         - None

        Stores:
         - "sigma_ug" in self

        Returns:
         - second moment of :math:`alpha`

        Mathematically performing the operation

        .. math:: \\frac{\\langle (\\ug - \\langle \\ug \\rangle)^2 \\rangle }{\\langle \\ug \\rangle^2}
        
        """

        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        ug_avg = self.area_avg('ug1')

        for angle, r_dict in self.data.items():
            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                    vars_temp.append(rstar * (midas_dict['ug1'] - ug_avg)**2)
                    if debug: print(angle, midas_dict, file=debugFID)
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            if debug: print("Arrays to integrate", rs, vars, file=debugFID)
                
            param_r.append( integrate.simpson(y=vars, x=rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(y=vars, x=rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)
        param_r_int = [var for _, var in sorted(zip(angles, param_r))]
        angles_int = sorted(angles)
        I = integrate.simpson(y=param_r_int, x=angles_int) / np.pi / ug_avg**2 # Integrate wrt theta, divide by normalized area
        if debug: print('Calculated sigma_ug: ', I)

        self.sigma_ug = I
        return I


    def calc_mu3_alpha(self):
        """Calculates the third moment of alpha

        Inputs:
         - None

        Stores:
         - "mu3_alpha" in self

        Returns:
         - Third moment of :math:`alpha`

        Mathematically performing the operation

        .. math:: \\frac{\\langle (\\alpha - \\langle \\alpha \\rangle)^3 \\rangle }{\\langle \\alpha \\rangle^3}
        
        """

        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        alpha_avg = self.area_avg('alpha')

        for angle, r_dict in self.data.items():
            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                    vars_temp.append(rstar * (midas_dict['alpha'] - alpha_avg)**3)
                    if debug: print(angle, midas_dict, file=debugFID)
            
            vars = [var for _, var in sorted(zip(rs_temp, vars_temp))]
            rs = sorted(rs_temp)

            if debug: print("Arrays to integrate", rs, vars, file=debugFID)
                
            param_r.append( integrate.simpson(y=vars, x=rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(y=vars, x=rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)
        param_r_int = [var for _, var in sorted(zip(angles, param_r))]
        angles_int = sorted(angles)
        I = integrate.simpson(y=param_r_int, x=angles_int) / np.pi / alpha_avg**3 # Integrate wrt theta, divide by normalized area
        if debug: print('Calculated mu3_alpha: ', I)

        self.mu3_alpha = I
        return I

    def top_bottom(self, param, even_opt='first') -> float:
        """Honestly, I forgot what this does

        I think it area-averages the way Bottin did, which is not
        a good way of doing so. 
        
        """
        
        # Check that the parameter that the user requested exists
        try:
            dummy = self.data[90][1.0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.data, file=debugFID)
            print(f"Could not area-average {param} for condition {self.name}")
            return
        
        # We have to integrate twice, once with resepect to r, again with respect to phi
        # Start with r

        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        self.mirror()

        for angle, r_dict in self.data.items():

            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0: # This should be unnecessary now with the new mirror, but it's not hurting anyone by being here
                    try:
                        rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                        vars_temp.append( midas_dict[param] * rstar)
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

        I = integrate.simpson(y=param_r, x=angles) / np.pi # Integrate wrt theta, divide by normalized area
        return I

    def calc_dpdz(self, method = 'LM', m = 0.316, n = 0.25, chisholm = 25, k_m = 0.10, L = None, alpha = None, akapower = 0.875):
        """ Calculates the frictional pressure gradient, dp/dz, according to various methods. Can access later with self.dpdz
 
        Options:
            - method    : Calculation method
              * 'LM'    : Lockhart Martinelli
              * 'Kim'   : Kim-modified Lockhart Martinelli
            - rho_f     : Liquid phase density
            - rho_g     : Gas phase density
            - mu_f      : Liquid phase dynamic viscosity
            - mu_g      : Gas phase dynamic viscosity
            - m         : Fed to calc_fric()
            - n         : Fed to calc_fric()
            - chisholm  : Chisholm parameter, the C in Lockhart-Martinelli
            - k_m       : Minor loss coefficient
            - L         : Length of restriction, only matters for 'Kim' method

        """

        if chisholm == 'Ryan':
            chisholm = 26 - 4.7*np.cos( self.theta*np.pi/180 )
        
        if method.lower() == 'lm' or method.lower() == 'lockhart' or method.lower() == 'lockhart-martinelli':
            f_f, f_g = self.calc_fric(m = m, n = n)

            dpdz_f = f_f * 1/self.Dh * self.rho_f * self.jf**2 / 2
            dpdz_g = f_g * 1/self.Dh * self.rho_g * self.jgloc**2 / 2
            chi2 = dpdz_f / dpdz_g

            phi_f2 = 1 + chisholm/np.sqrt(chi2) + 1 / chi2
            dpdz = phi_f2 * dpdz_f

        elif method.lower() == 'kim':
            f_f, f_g = self.calc_fric(m = m, n = n)

            dpdz_f = f_f * 1/self.Dh * self.rho_f * self.jf**2 / 2
            dpdz_g = f_g * 1/self.Dh * self.rho_g * self.jgloc**2 / 2
            dpdz_m = k_m * self.rho_f * self.jf**2 / 2 / L

            chi2 = dpdz_f / dpdz_g
            chiM2 = dpdz_f / dpdz_m

            phi_f2 = (1 + 1 / chiM2) + np.sqrt(1 + 1 / chiM2) * chisholm / np.sqrt(chi2) + 1 / chi2
            dpdz = phi_f2 * dpdz_f

        elif method.lower() == 'akagawa':
            f_f, f_g = self.calc_fric(m = m, n = n)

            dpdz_f = f_f * 1/self.Dh * self.rho_f * self.jf**2 / 2
            phi_f2 = ((1 - alpha)**(-akapower))**2

            dpdz = phi_f2 * dpdz_f

        else:
            raise NotImplementedError(f'{method} is not a valid option for calc_dpdz. Try "LM" ')
        
        self.dpdz = dpdz
        self.tau_w = dpdz * self.Dh / 4

        return dpdz

    def calc_vwvg(self) -> None:
        """Calculates void weighted Vgj

        uses self.jgloc / self.area_avg('alpha'), which is not a great method

        """

        print("This guy needs work, probably don't want to use it")
        self.vwvg = self.jgloc / self.area_avg('alpha')
        return self.jgloc / self.area_avg('alpha')
    
    def calc_W(self):
        """Calculates the third moment of alpha

        Inputs:
         - None

        Stores:
         - "W" in midas_dict

        Returns:
         - Area-average W

        Mathematically performing the operation

        .. math:: W = \\frac{v_r}{v_f}
        
        """

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    vfp = midas_dict['vf']
                except KeyError:
                    self.calc_vr()

                if vfp == 0:
                    if midas_dict['vr'] != 0:
                        print("Warning, vf = 0 but vr != 0. Defaulting to W = 0")
                    midas_dict.update({
                            'W' : 0
                        })
                else:

                    midas_dict.update({
                        'W' : midas_dict['vr'] / vfp
                    })

        return self.area_avg('W')
    
    def calc_fric(self, method = 'Blasius', m = 0.316, n=0.25):
        """ Calculates friction factor for each phase based on bulk Re

        .. math:: Re_k = \\frac{\\rho_k j_k D}{\\mu_k}

        Method Options:
         - Blasius
        .. math:: f_k = \\frac{C}{Re_k^n}

        Returns:
         - Tuple of f_f, f_g
        
        Stores
         -self.ff
         -self.fg
         -self.tau_fw
        
        """
        
        Re_f = self.rho_f * self.jf * self.Dh / self.mu_f
        Re_g = self.rho_g * self.jgloc * self.Dh / self.mu_g

        if method.lower() == 'blasius':
            ff = m / Re_f**n
            fg = m / Re_g**n
        else:
            raise NotImplementedError("Invalid method, try Blasius")

        self.ff = ff
        self.fg = fg
        tau_fw = self.ff/4 * self.rho_f * self.jf**2/2
        self.tau_fw = tau_fw
        return (ff, fg)

    
    def calc_mu_eff(self, method='Ishii', alpha_peak = 1.0):
        """Method for effective/mixture viscosity
        
        Inputs:
         - method, what method to use for modeling :math:`\\mu_{eff}`
         - alpha_peak, parameter for some models
        
        Stores:
         - "mu_eff" in midas_dict 
         - "mu_m" in midas_dict

        Options for method:
         - Ishii, :math:`\\mu_{eff} = \\mu_{m} = \\mu_{f} (1 - \\frac{\\alpha}{ \\alpha_{max} })^{-2.5 \\alpha_{max} \\frac{\\mu_g + 0.4 \\mu_g}{\\mu_g + \\mu_g} }`
         - Ishii_AA, :math:`\\mu_{eff} = \\mu_{m} = \\mu_{f} (1 - \\frac{\\langle \\alpha \\rangle}{ \\alpha_{max} })^{-2.5 \\alpha_{max} \\frac{\\mu_g + 0.4 \\mu_g}{\\mu_g + \\mu_g} }`

        Returns:
         - area average effective viscosity

        """

        self.mirror()
        alpha_avg = self.area_avg('alpha')
                
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():

                if method.lower() == 'ishii':
                    mu_m = self.mu_f * (1 - midas_dict['alpha'] / alpha_peak)**(-2.5*alpha_peak * (self.mu_g + 0.4*self.mu_f) / (self.mu_g + self.mu_f)  )

                elif method.lower() == 'ishii_AA':
                    mu_m = self.mu_f * (1 - alpha_avg / alpha_peak)**(-2.5*alpha_peak * (self.mu_g + 0.4*self.mu_f) / (self.mu_g + self.mu_f)  )
                    
                elif method.lower() == 'avg_void':
                    mu_m = self.mu_f / (1 - alpha_avg)

                mu_eff = mu_m
                mu_m = mu_eff

                midas_dict.update({'mu_eff': mu_eff})
                midas_dict.update({'mu_m': mu_m})

        return self.area_avg('mu_eff')

    def calc_cd(self, method='Ishii-Zuber', vr_cheat = False, limit = 0, const_CD = 0.44):
        """Method for calculating drag coefficient
        
        Inputs:
         - method, what method to use for modeling :math:`C_{D}`
         - vr_cheat, flag to use "vr" from midas_dict or "vr_model" when calculating :math:`Re_{b}`
         - limit, if supplied, will limit the drag coefficient to the given maximum value. For instance, 0.44 like in CFX
        
        Stores:
         - "cd" in midas_dict 
         - "Reb" in midas_dict. Calculated by :math:`Re_{b} = \\frac{(1 - \\alpha) \\rho_{f} v_{r} D_{sm,1} }{\\mu_{m}}`. \ 
:math:`\\mu_{m}` comes from :any:`calc_mu_eff`
         - "eo", if using "tomiyama" or "ishii" limits

        Options for method:
         - Ishii-Zuber, :math:`C_{D} = \\frac{24}{Re_{b}} (1 + 0.1 Re_{b}^{0.75})`
         - Schiller-Naumann, :math:`C_{D} = \\frac{24}{Re_{b}} (1 + 0.15 Re_{b}^{0.687})`. Here, :math:`Re_{b}` uses :math:`\\mu_{f}`
         - const

        Returns:
         - area average drag coefficient

        """

        self.calc_mu_eff()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():

                if vr_cheat:
                    Reb = (1 - midas_dict['alpha']) * midas_dict['Dsm1'] * self.rho_f * abs(midas_dict['vr']) / midas_dict['mu_m']
                else:

                    if 'vr_model' not in midas_dict.keys(): # Initialize for iteration
                        midas_dict.update({'vr_model': -1})

                    Reb = (1 - midas_dict['alpha']) * midas_dict['Dsm1'] * self.rho_f * abs(midas_dict['vr_model']) / midas_dict['mu_m']

                midas_dict.update({'Reb': Reb})

                if method == 'Ishii-Zuber' or method == 'IZ' or method == 'Ishii':

                    if Reb > 0:
                        cd = 24/Reb * (1 + 0.1*Reb**0.75)
                    else:
                        cd = 0

                elif method.lower() == 'schiller-naumann':
                    Reb = (1 - midas_dict['alpha']) * midas_dict['Dsm1'] * self.rho_f * abs(midas_dict['vr_model']) / self.mu_f

                    cd = 24/Reb * (1 + 0.15*Reb**0.687)
                
                elif method == 'constant' or method == 'const':
                    cd = const_CD

                if type(limit) == str:
                    if limit.lower() == "tomiyama":
                        eo = self.g * (self.rho_f - self.rho_g) * midas_dict['Dsm2']
                        limit = 8/3 * eo / (eo + 4)
                        midas_dict.update({'eo': eo})
                        
                    elif limit.lower() == 'ishii-chawla':
                        eo = self.g * (self.rho_f - self.rho_g) * midas_dict['Dsm2']
                        limit = min(2/3*np.sqrt(eo), 8/3)
                        midas_dict.update({'eo': eo})
                    
                    else:
                        raise NotImplementedError(f"{limit} not a valid type for limiting behavior. Please enter tomiyama, Ishii-Chawla, or set a constant limit (e.g. limit = 0.44)")
                
                cd = max(limit, cd) # Either 0, set by user, or set by above string

                midas_dict.update({'cd': cd})

        return self.area_avg('cd')

    def calc_vr_model(self, method='km1_simp', kw = -0.98, n=1, Lw = 5, kf = 0.089, 
                      iterate_cd = True, initial_vr = None, 
                      quiet = True, recalc_cd = True, iter_tol = 1e-4, custom_f = None, CC = 1):
        """Method for calculating relative velocity based on models
        
        Inputs:
         - method, what method to use for modeling :math:`v_{r}`
         - kw, wake coeffieient for some models
         - n, exponent for some models
         - Lw, really :math:`L_{w}^{*}`, effective wake length divided by the bubble diameter
         - kf, fluid coefficient for some models
         - iterate_cd, flag to pass to :any:`calc_cd`. Basically whether to iterate to calculate :math:`C_{D}` based on \
the newly calculated :math:`v_{r}` or not
         - quiet, flag for extra debugging messages
        
        Stores:
         - "vr\_'method_name'" in midas_dict 
         - "vr_model" in midas_dict 

        Options for method:
         - "km1_simp", most up to date :math:`v_{r} = K_{f} - K_{w}  \\alpha  v_{f}  C_{d}^{1/3}`
         - "Ishii-Chawla", :math:`v_{r} = \\sqrt{2} (\\frac{\\sigma g \\Delta \\rho}{\\rho_{f}^{2}})^{1/4}`
         - A bunch of obsolete ones that I haven't deleted but probably should

        Returns:
         - area average relative velocity calculated by the model

        """

        MAX_ITERATIONS = 1000
        iterations = 0
        initialize_vr = True

        if initial_vr is None:
            if self.theta == 90:
                initial_vr = 0.5
            else:
                initial_vr = -0.5

        while True:
            if recalc_cd:
                if iterate_cd:
                    
                    if initialize_vr:
                        for angle, r_dict in self.data.items():
                            for rstar, midas_dict in r_dict.items():
                                midas_dict.update(
                                    {'vr_model': initial_vr}
                                )
                        initialize_vr = False

                    self.calc_cd(vr_cheat=False)
                else:
                    self.calc_cd(vr_cheat=True)
            try:
                old_vr = self.area_avg('vr_model', recalc=True)
            except KeyError:
                old_vr = initial_vr # Initialize?

            vr_name = "vr_" + method

            for angle, r_dict in self.data.items():
                for rstar, midas_dict in r_dict.items():
                    
                    if method == 'wake_1':
                        vr = kw  * midas_dict['vf'] * midas_dict['cd']**(1./3)
                    
                    elif method == 'wake_alpha':
                        vr = kw  * (1 - midas_dict['alpha'])**n * midas_dict['vf'] * midas_dict['cd']**(1./3)
                    
                    elif method == 'wake_alpha2':
                        vr = kw  * (midas_dict['alpha']*(1 - midas_dict['alpha']))**n * midas_dict['vf'] * midas_dict['cd']**(1./3)

                    elif method == 'wake_lambda':
                        self.calc_avg_lat_sep()
                        vr = kw  * midas_dict['vf'] * midas_dict['cd']**(1./3) * midas_dict['Dsm1']**(2/3) * midas_dict['lambda']**(-2./3)

                    elif method == 'wake_vg_lambda':
                        self.calc_avg_lat_sep()
                        vr = midas_dict['ug1'] /(1 + kw * midas_dict['cd']**(1./3) * midas_dict['Dsm1']**(2/3) * midas_dict['lambda']**(-2./3) )

                    elif method == 'hubris' or method == 'wake_vr':
                        cd = midas_dict['cd']
                        c2 = kw
                        vr = 1/2* midas_dict['ug1']*(1 - 3 *kw* cd**(1/3)*(2**(2/3) - 2* Lw**(1/3)) - 2* Lw + 6* c2**(3/2)* np.sqrt(cd) *(np.arctan(1./(2**(1/3) * np.sqrt(c2) * cd**(1/6))) - np.arctan(1/((np.sqrt(c2)* cd**(1/6))/Lw**(1/3)))))
                    
                    elif method == 'km1':
                        vr = kw * (np.pi/4)**(1/3) * midas_dict['alpha'] * midas_dict['vf'] * midas_dict['cd']**(1./3) *  (2**(-1./3) - Lw**(1/3))/(0.5 - Lw) + kf * midas_dict['vf']
                    
                    elif method == 'km1_simp' or method == 'prelim':
                        vr = -kw * midas_dict['alpha'] * midas_dict['vf'] * midas_dict['cd']**(1./3) - kf * midas_dict['vf']

                    elif method == 'prelim_plus':
                        if custom_f is None:
                            ff, fg = self.calc_fric()
                        else:
                            ff = custom_f

                        # if midas_dict['Dsm1'] == 0 and rstar != 1.0:
                        #     print(self, angle, rstar)

                        try:
                            vr = (
                            -kw * midas_dict['alpha'] * midas_dict['vf'] * midas_dict['cd']**(1./3) - kf * midas_dict['vf'] 
                            + np.sqrt( 4./3 * self.void_area_avg('Dsm1')*0.001/midas_dict['cd'] * ( ff/self.Dh * self.jf**2/2 + 
                                                                                 (1 - midas_dict['alpha'])*(1-self.rho_g/self.rho_f) * self.gz ) )
                            )
                        except ZeroDivisionError:
                            vr = 0

                        if vr == np.inf and midas_dict['cd'] == 0:
                            if not quiet: warnings.warn(f"vr nan for {angle, rstar}, setting to 0")
                            vr = 0
                        
                        if vr != vr:
                            if not quiet: warnings.warn(f"vr nan for {angle, rstar}, setting to 0")
                            vr = 0
                        

                    elif method == 'proper_integral':
                        warnings.warn("This method is probably no good, messed up the math")
                        vr = midas_dict['ug1'] / ( (0.5 - Lw) + kw * midas_dict['cd']**(1./3) *(np.pi/4)**(1/3)* (2**(-1./3) - Lw**(1/3)))

                    elif method == 'proper_integral_alpha':
                        warnings.warn("This method is probably no good, messed up the math")
                        vr = midas_dict['ug1'] / ( (0.5 - Lw) + kw * midas_dict['alpha']**n *midas_dict['cd']**(1./3) *(np.pi/4)**(1/3)* (0.5**(1./3) - Lw**(1/3)))

                    elif method.lower() == 'ishii-chawla' or method.lower() == 'ishii' or method.lower() == 'ishii chawla':
                        vr = np.sqrt(2) * (self.sigma * self.g * (self.rho_f - self.rho_g) / (self.rho_f**2))**0.25
                    
                    elif method.lower() == 'chahed':
                        self.calc_dpdz()
                        self.calc_grad('alpha')
                        self.calc_grad('vf')
                        if abs(midas_dict['alpha']) < 1e-6 or abs(midas_dict['cd']) < 1e-6:
                            discrim = 0
                        else:
                            # print(self.Dh, self.rho_f, midas_dict['cd'], midas_dict['alpha'])
                            discrim = 4/3 * midas_dict['Dsm1'] / midas_dict['cd'] * ( 4/self.Dh * self.tau_w  / self.rho_f + self.gz*(1-midas_dict['alpha']) - CC/midas_dict['alpha'] * midas_dict['grad_alpha_r'] * midas_dict['grad_vf_r'])

                        vr = float(np.sign(discrim) * np.sqrt(discrim))
                    else:
                        print(f"{method} not implemented")
                        return -1
            
                    midas_dict[vr_name] = vr
                    midas_dict['vr_model'] = vr


            iterations += 1

            try:
                self.area_avg('vr_model', recalc=True)
            except:
                print(f"Error calculating vr" )

            if old_vr == np.inf or old_vr != old_vr or vr == np.inf:
                print(f"Error calculating vr: {self.area_avg('vr_model', recalc=True)}")

            if old_vr == 0:
                if not quiet:
                    print(f"vr_model calculated as 0 after {iterations} iterations")
                    print(old_vr, self.area_avg('vr_model', recalc=True))
                return

            if abs(old_vr - self.area_avg('vr_model', recalc=True)) / abs(old_vr) < iter_tol:
                if not quiet:
                    print(f"vr_model converged in {iterations} iterations")
                    print(old_vr, self.area_avg('vr_model', recalc=True))
                return self.area_avg("vr_model")
            
            if iterations > MAX_ITERATIONS:
                print("Warning, max iterations exceeded in calculating vr_model")
                print(f"{old_vr}\t{self.area_avg('vr_model', recalc=True)}\t{(old_vr - self.area_avg('vr_model', recalc=True))/old_vr*100}")
                return
            
        
        return self.area_avg("vr_model")

            
    def calc_vgj_model(self):
        """Method for calculating local Vgj based on models
        
        midas_dict['vgj_model'] = (1 - midas_dict['alpha']) * midas_dict['vr_model']

        """

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():

                if 'vr_model' not in midas_dict.keys():
                    print(f"Warning: vr_model not found for {angle}, {rstar}, calling calc_vr_model with default inputs")
                    self.calc_vr_model()

                midas_dict['vgj_model'] = (1 - midas_dict['alpha']) * midas_dict['vr_model']

        return
    
    def calc_IS_term(self, method = 'power', n=2, mu = 1.5):
        self.calc_cd()
        self.calc_fric()
        self.calc_grad('alpha')
        
        for angle, r_dict in self.data.items():
            for r_star, midas_dict in r_dict.items():
                    
                    if r_star == 0 and n < 0:
                        midas_dict.update(
                            {'IS': 0 , 'ISxgrad': 0}
                         )
                        continue

                    if method == 'power':
                        taui = float(self.tau_fw * r_star**n)
                        midas_dict.update(
                            {'ISxgrad': float(taui * ( midas_dict['grad_alpha_r'] )) , 'IS': taui}
                         )
        
                    elif method == 'power_total':
                        taui = float(self.tau_fw * r_star**n)
                        midas_dict.update(
                            {'ISxgrad': float(taui * (midas_dict['grad_alpha_total'] )) , 'IS': taui}
                        )
                        
                    elif method == 'lognorm':
                        taui = self.tau_fw * np.exp( -(np.log(1-r_star) + mu)**2 )
                        midas_dict.update(
                            {'ISxgrad': float(taui * (midas_dict['grad_alpha_r'] )) , 'IS': float(taui)}
                        )

                    elif method == 'alpha':
                        taui = n * self.tau_fw * midas_dict['alpha'] / self.area_avg('alpha')
                        midas_dict.update(
                            {'ISxgrad': float(taui * (midas_dict['grad_alpha_r'] )) , 'IS': float(taui)}
                        )


        return self.area_avg('ISxgrad')
    
    def calc_aa_vr_model(self, method='km1_naive', IS_method = 'power', kw=-0.98, kf=-0.083, Lw = 5, Cavf=1, Ctau=1, n=2, IS_mu = 1.5, Cvracd = 1):

        if method == 'km1_naive':
            vr = kw * (np.pi/4)**(1/3) * self.area_avg('alpha') * self.jf / (1 - self.area_avg('alpha')) * self.area_avg('cd')**(1./3) *  (2**(-1./3) - Lw**(1/3))/(0.5 - Lw) + kf * self.jf / (1 - self.area_avg('alpha'))
        
        elif method == 'km1_naive2' :#or method == 'prelim':
            self.calc_cd()
            vr = kw * self.area_avg('alpha') * self.jf / (1 - self.area_avg('alpha')) * self.area_avg('cd')**(1./3)  + kf * self.jf / (1 - self.area_avg('alpha'))

        elif method == 'km1_simp' or method == 'prelim':
            vr = -kw * self.area_avg('alpha') * self.jf / (1 - self.area_avg('alpha')) * self.area_avg('cd')**(1./3) - kf * self.jf / (1 - self.area_avg('alpha'))
        
        elif method == 'final':
            vr = Cvracd * -kw * self.area_avg('alpha') * self.jf / (1 - self.area_avg('alpha')) * self.area_avg('cd')**(1./3) - kf * self.jf / (1 - self.area_avg('alpha'))
        
        elif method == 'IS_Ctau':
            self.calc_cd()
            self.calc_fric()
            rb = self.void_area_avg('Dsm1') / 2 /1000 # Convert to m
            CD = self.void_area_avg('cd')
            alpha = self.area_avg('alpha')

            discrim = (1-alpha)*self.gz * (self.rho_f - self.rho_g) + (1-Ctau)*4*self.tau_fw/self.Dh
            vr = np.sign(discrim) * np.sqrt(8*rb/3 * 1/(CD * self.rho_f) * abs( discrim ))

        elif method == 'IS':
            IS_term = self.calc_IS_term(method = IS_method, n = n, mu = IS_mu)

            rb = self.void_area_avg('Dsm1') / 2 /1000 # Convert to m
            CD = self.void_area_avg('cd')
            alpha = self.area_avg('alpha')

            discrim = 4*self.tau_fw/self.Dh + (1-alpha)*self.gz * (self.rho_f - self.rho_g) - 1/alpha * IS_term
            vr = np.sign(discrim) * np.sqrt(8*rb/3 * 1/(CD * self.rho_f) * abs( discrim ))
        else:
            raise(ValueError("Invalid calc_aa_vr_model method"))
            

        self.vwvgj = (1-self.area_avg('alpha'))*vr
        self.aa_vr = vr
        return vr
    
    def calc_errors(self, param1:str, param2:str) -> float:
        """ Calculates the errors, ε, between two parameters (param1 - param2) in midas_dict

        Usually want to do param1=predicted, param2=experimental
        
        Stores:
         - error, "eps_param1_param2", param1 - param2
         - relative, "eps_rel_param1_param2", (param1 - param2) / param2
         - absolute relative, "eps_abs_rel_param1_param2", | param1 - param2 | / param2
         - square, "eps_sq_param1_param2", (param1 - param2)**2
         - relative square, "eps_rel_sq_param1_param2", ((param1 - param2)/param2)**2

         If param2 = 0, relative errors are considered 0
        
        """

        param_error_name = "eps_" + param1 + "_" + param2
        param_rel_error_name = "eps_rel_" + param1 + "_" + param2
        param_abs_rel_error_name = "eps_abs_rel_" + param1 + "_" + param2
        param_sq_error_name = "eps_sq_" + param1 + "_" + param2
        param_rel_sq_error_name = "eps_rel_sq_" + param1 + "_" + param2

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                midas_dict[param_error_name] = midas_dict[param1] - midas_dict[param2]
                midas_dict[param_sq_error_name] = (midas_dict[param1] - midas_dict[param2])**2

                if midas_dict[param2] != 0:
                    midas_dict[param_rel_error_name] = (midas_dict[param1] - midas_dict[param2]) / midas_dict[param2]
                    midas_dict[param_rel_sq_error_name] = ((midas_dict[param1] - midas_dict[param2]) / midas_dict[param2])**2
                    midas_dict[param_abs_rel_error_name] = abs(midas_dict[param1] - midas_dict[param2]) / abs(midas_dict[param2])
                else:
                    midas_dict[param_rel_error_name] = 0
                    midas_dict[param_rel_sq_error_name] = 0
                    midas_dict[param_abs_rel_error_name] = 0

        return self.area_avg(param_error_name)
    
    def calc_AA_error(self, param1:str, param2:str) -> float:
        """ Calculates the error, ε, between the area-average of two parameters (⟨param1⟩ - ⟨param2⟩) in midas_dict

        Usually want to do param1=predicted, param2=experimental
        
        Returns
         - relative error (⟨param1⟩ - ⟨param2⟩) / ⟨param1⟩
        
        """
        eps = self.area_avg(param1) - self.area_avg(param2)
        rel_error = eps / self.area_avg(param1)

        return rel_error

    def calc_symmetry(self, param, sym_type = 'sym90', method = 'rmse', rel_error = False):
        """ Function for checking the symmetry of a condition

        sym_type
         - 'sym90', will compare data at 0° and 180°

        .. math:: \\epsilon(r) = param(r, 0°) - param(r, 180°)

        rel_error will divide by param 0

        Method
         - rmse, 
        
        .. math:: \\epsilon_{sym} = \\sqrt{\\bar{\\epsilon(r)^{2}}}

        Returns
         - :math:`\\epsilon_{sym}`
        
        Stores
         - :math:`\\epsilon_{sym}`
        
        """
        errs = []
        if sym_type == 'sym90':

            for rstar, r_dict in self.data[0].items():
                try:
                    err = r_dict[param] - self.data[180][rstar][param]
                except KeyError:
                    continue
                if rel_error:
                    try:
                        err = err / r_dict[param]
                    except ZeroDivisionError:
                        err = 0
                errs.append(err)

        errs = np.asarray(errs)

        if method == 'rmse':
            sym_error = np.sqrt(np.mean(errs**2))
        elif method == 'mean':
            sym_error = np.mean(errs)

        self.sym_error = sym_error

        return self.sym_error

    def calc_symmetry_area_avg(self, param, sym_type='sym_half', rel_error=True, even_opt='first'):
        """ Function for checking the area-averaged symmetry of a condition between the left and right half.
    
        Inputs:
        - param: The parameter for which to calculate the symmetry
        - sym_type: 'sym_half' to compare left half (0-180) and right half (180-360)
        - rel_error: If true, the error will be relative to the parameter at the left half angles
        - even_opt: Option for integration method in simpson's rule (as used in area_avg)
        Returns:
       - Area-averaged symmetry error across all radial and angular positions
       """
        
        if sym_type != 'sym_half':
           raise ValueError("This function is only designed for 'sym_half' symmetry type.")

         # Define angle pairs for left-right half symmetry
        angle_pairs = [(0, 180), (22.5, 202.5), (45, 225), (67.5, 247.5), (90, 270), (292.5, 112.5), (315, 135), (337.5, 157.5)]

        sym_errors = []  # To store symmetry errors for each angle pair
        areas = []       # To store areas for weighted integration

         # Iterate over angle pairs
        for angle1, angle2 in angle_pairs:
            if angle1 not in self.data or angle2 not in self.data:
                continue  # Skip if angles are not in the data

            param_r = []  # To store the radial integration result
            rs_list = []  # Radial positions to integrate over
            
            # Calculate symmetry error for each radial position
            for rstar, r_dict in self.data[angle1].items():
                try:
                  sym_error = abs(r_dict[param] - self.data[angle2][rstar][param])
                except KeyError:
                   continue  # Skip if data is missing for the corresponding rstar in angle2

                if rel_error:
                    try:
                       sym_error /= r_dict[param]  # Calculate relative error
                    except ZeroDivisionError:
                       sym_error = 0  # Handle division by zero case
            
                rs_list.append(rstar)
                param_r.append(abs(sym_error) * rstar)  # Multiply by rstar for radial area weighting

            # Sort by radial position for integration
            param_r_sorted = [x for _, x in sorted(zip(rs_list, param_r))]
            rs_sorted = sorted(rs_list)
            
            # Integrate wrt r (radial direction) using simpson's rule
            if len(rs_sorted) > 1:
                try:
                    radial_avg_sym = integrate.simpson(y=param_r_sorted, x=rs_sorted)
                except Exception as e:
                    print(f"Integration error: {e}")
                    radial_avg_sym = 0
            else:
                radial_avg_sym = 0

            sym_errors.append(radial_avg_sym)
            areas.append(1)  # For now, give equal weighting to all angles

            # Now integrate over angles (θ)
        angles_in_radians = [angle_pair[0] * np.pi / 180 for angle_pair in angle_pairs]
        if len(sym_errors) > 1:
            try:
               # Final integration wrt angle θ
              area_avg_sym_error = integrate.simpson(y=sym_errors, x=angles_in_radians) / np.pi
            except Exception as e:
                print(f"Angle integration error: {e}")
                area_avg_sym_error = 0
        else:
           area_avg_sym_error = 0

        return area_avg_sym_error


    def calc_avg_lat_sep(self):
        """Calculates average lateral separation distance between bubbles

        Inputs:
         - None

        Stores:
         - "lambda" in midas_dict
         - "lambda*" in midas_dict
         - "alpha_lambda" in midas_dict

        Returns:
         - Area-average :math:`\\lambda`

        Assuming

        .. math:: \\lambda = \\frac{v_g}{f}

        Where :math:`f` is the bubble frequency. Also estimates the void fraction based on this :math:`\\lambda`,

        .. math:: \\alpha = \\frac{D_{sm,1}}{\\lambda}
        
        """

        self.mirror()

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                
                try:
                    midas_dict['lambda'] = midas_dict['ug1'] / midas_dict['bub_freq']
                    midas_dict['lambda*'] = midas_dict['lambda'] / (midas_dict['Dsm1']/1000)
                    midas_dict['alpha_lambda'] = midas_dict['Dsm1']/1000 / midas_dict['lambda']
                except ZeroDivisionError:
                    midas_dict['lambda'] = np.inf
                    midas_dict['lambda*'] = np.inf
                    midas_dict['alpha_lambda'] = 0

        return self.area_avg('lambda')
    

    def calc_COV_RC(self, alpha_peak = 0.75, alpha_cr = 1.00, avg_method = 'legacy', reconstruct_flag = False, debug = False):
        """Calculates the experimental Random Collision Covariance based on Talley (2012) method (without modification factor m_RC)
         - Stored in self.COV_RC
         
         Inputs:
         - alpha_peak, maximum void fraction based on hexagonal-closed-packed (HCP) bubble distribution
         - alpha_cr, critical alpha to activate Random Collision, Talley (2012), Kong (2018) 
        
         Authors:
         - Quan (05/15/2024)
         - Kang (08/19/2024)
        """

        if debug:
            print(f"\t_________________________________________________________")
            print(f"\t                         COV_RC                          ")
            print(f"\tjf = {self.jf}, jgref = {self.jgref}, L/D = {self.LoverD}")

        # Identical to calc_COV_TI
        ########################################################################################################################
        if reconstruct_flag == True:
            alpha_str = 'alpha_reconstructed'
        else:
            alpha_str = 'alpha'

        # Temporary, replace later
        alpha_avg       = self.area_avg('alpha',method=avg_method)
        alpha_avg       = round(alpha_avg,3)

        ai_avg          = self.area_avg('ai',method=avg_method)
        ai_avg          = round(ai_avg,2)

        # rho_f           = self.rho_f                                        # Liquid phase density [kg/m**3]
        # rho_g           = self.rho_g                                        # Gas phase density [kg/m**3]
        # mu_f            = self.mu_f                                         # Dynamic viscosity of water [Pa-s]
        # sigma           = self.sigma                                        # Surface tension

        rho_f           = 998                                                # What Quan used
        rho_g           = 1.226                                              # What Quan used
        mu_f            = 0.001
        sigma           = self.sigma
        Dh              = self.Dh                                           # Hydraulic diameter [m]
                
        rho_m           = (1 - alpha_avg) * rho_f + alpha_avg * rho_g       # Mixture density
        mu_m            = mu_f / (1 - alpha_avg)                            # Mixture viscosity
        v_m             = (rho_f * self.jf + rho_g * self.jgloc) / rho_m    # Mixture velocity
        # Rem             = rho_m * v_m * Dh / mu_m                         # Mixture Reynolds number, ***CAREFUL*** I have seen some versions of the IATE script that use rho_f instead as an approximation
        Rem             = rho_m * v_m * Dh / mu_f                           # Mixture Reynolds number, older versions use mu_f as an approximation
        f_TP            = 0.316 * (mu_m / mu_f / Rem)**0.25                 # Two-phase friction factor, Talley (2012) and Worosz (2015), also used in iate_1d_1g
        eps             = f_TP * v_m**3 / 2 / Dh                            # Energy dissipation rate (Wu et al., 1998; Kim, 1999), also used in iate_1d_1g        
        eps             = round(eps,2)

        # Switch away from using data
        alpha_avg       = self.area_avg(alpha_str,method=avg_method)

        Dsm_exp         = 1000 * 6 * alpha_avg / ai_avg             # for void re-construction method only 
        Dsm_exp         = round(Dsm_exp,2) / 1000

        if debug:
            print(f"\t\trho_m: {rho_m}\tmu_m: {mu_m}\tv_m: {v_m}\tRem: {Rem}\tf_TP: {f_TP}\teps: {eps}")
            print(f"\t\tDsm_exp: {Dsm_exp}\talpha_str_avg: {alpha_avg}\tai_avg: {ai_avg}")

        for angle, r_dict in self.data.items():
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
                    u_t_loc = 1.4 * np.cbrt(eps) * np.cbrt(Db_loc)              # Turbulent velocity (Batchelor, 1951; Rotta, 1972), also used in iate_1d_1g
                else:
                    u_t_loc = 0                                                 # TI and RC are driven by the turbulent fluctuation velocity (u_t)

                midas_dict['u_t_loc'] = u_t_loc     #Quan, 1106

                ########################################################################################################################
                # Talley 2012, section 3.3.1
                COV_RC_loc = u_t_loc * ai_loc**2 / (np.cbrt(alpha_peak) * (np.cbrt(alpha_peak) - np.cbrt(alpha_loc)))
                
                midas_dict['COV_RC_loc'] = COV_RC_loc

                if debug:
                    # print(f"\t\t\t{angle:2.1f}\t{rstar:.2f}\t|\talpha: {alpha_loc:.4f}\tDbloc: {Db_loc:.4f}\tCOV_RC_loc: {COV_RC_loc:.4f}")
                    pass
                
        # Talley does not area-average local u_t; instead computes <u_t> with area-averaged parameters
        u_t_avg = 1.4 * np.cbrt(eps) * np.cbrt(6 * alpha_avg / ai_avg)       
       # u_t_avg = self.area_avg('u_t_loc', method=avg_method)   #Quan, 1106

        if u_t_avg > 0:
            COV_RC_avg = u_t_avg * ai_avg**2 / (np.cbrt(alpha_peak) * (np.cbrt(alpha_peak) - np.cbrt(alpha_avg)))
            COV_RC = self.area_avg('COV_RC_loc',method=avg_method) / COV_RC_avg

            if debug:
                print(f"\n\t\tu_t_avg: {u_t_avg}")
                print(f"\n\tCOV_RC: {COV_RC}\tCOV_RC_num: {self.area_avg('COV_RC_loc',method=avg_method)}\tCOV_RC_den: {COV_RC_avg}")
        else:
            COV_RC = 0

        self.COV_RC = COV_RC

        return COV_RC


    def calc_COV_TI(self, alpha_peak = 0.75, alpha_cr = 1.00, We_cr = 6, avg_method = 'legacy', reconstruct_flag = False, debug = False):

        if debug:
            print(f"\t_________________________________________________________")
            print(f"\t                         COV_TI                          ")
            print(f"\tjf = {self.jf}, jgref = {self.jgref}, L/D = {self.LoverD}")

        # Identical to calc_COV_RC
        ########################################################################################################################
        if reconstruct_flag == True:
            alpha_str = 'alpha_reconstructed'
        else:
            alpha_str = 'alpha'

        # Temporary, replace later
        alpha_avg       = self.area_avg('alpha',method=avg_method)
        alpha_avg       = round(alpha_avg,3)

        ai_avg          = self.area_avg('ai',method=avg_method)
        ai_avg          = round(ai_avg,2)

        # rho_f           = self.rho_f                                        # Liquid phase density [kg/m**3]
        # rho_g           = self.rho_g                                        # Gas phase density [kg/m**3]
        # mu_f            = self.mu_f                                         # Dynamic viscosity of water [Pa-s]
        # sigma           = self.sigma                                        # Surface tension

        rho_f           = 998                                                # What Quan used
        rho_g           = 1.226                                              # What Quan used
        mu_f            = 0.001
        sigma           = self.sigma
        Dh              = self.Dh                                           # Hydraulic diameter [m]
                
        rho_m           = (1 - alpha_avg) * rho_f + alpha_avg * rho_g       # Mixture density
        mu_m            = mu_f / (1 - alpha_avg)                            # Mixture viscosity
        v_m             = (rho_f * self.jf + rho_g * self.jgloc) / rho_m    # Mixture velocity
        # Rem             = rho_m * v_m * Dh / mu_m                           # Mixture Reynolds number, ***CAREFUL*** I have seen some versions of the IATE script that use rho_f instead as an approximation
        Rem             = rho_m * v_m * Dh / mu_f                           # Mixture Reynolds number, older versions use mu_f as an approximation
        f_TP            = 0.316 * (mu_m / mu_f / Rem)**0.25                 # Two-phase friction factor, Talley (2012) and Worosz (2015), also used in iate_1d_1g
        eps             = f_TP * v_m**3 / 2 / Dh                            # Energy dissipation rate (Wu et al., 1998; Kim, 1999), also used in iate_1d_1g        
        eps             = round(eps,2)

        # Switch away from using data
        alpha_avg       = self.area_avg(alpha_str,method=avg_method)

        Dsm_exp         = 1000 * 6 * alpha_avg / ai_avg
        Dsm_exp         = round(Dsm_exp,2) / 1000

        if debug:
            print(f"\t\trho_m: {rho_m}\tmu_m: {mu_m}\tv_m: {v_m}\tRem: {Rem}\tf_TP: {f_TP}\teps: {eps}")
            print(f"\t\tDsm_exp: {Dsm_exp}\talpha_str_avg: {alpha_avg}\tai_avg: {ai_avg}")

        for angle, r_dict in self.data.items():
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
                    u_t_loc = 1.4 * np.cbrt(eps) * np.cbrt(Db_loc)              # Turbulent velocity (Batchelor, 1951; Rotta, 1972), also used in iate_1d_1g
                else:
                    u_t_loc = 0                                                 # TI and RC are driven by the turbulent fluctuation velocity (u_t)

                midas_dict['u_t_loc'] = u_t_loc     #Quan, 1106

                ########################################################################################################################
                We = rho_f * u_t_loc**2 * Db_loc / sigma                        # Weber number criterion

                # Talley 2012, section 3.3.1
                if We >= We_cr:
                    COV_TI_loc = (u_t_loc * ai_loc**2 / alpha_loc) * np.sqrt(1 - (We_cr / We)) * np.exp(-We_cr / We)
                else:
                    COV_TI_loc = 0

                midas_dict['COV_TI_loc'] = COV_TI_loc

                if debug:
                    print(f"\t\t\t{angle:2.1f}\t{rstar:.2f}\t|\talpha: {alpha_loc:.4f}\tCOV_TI_loc: {COV_TI_loc:.4f}\tu_t: {u_t_loc}\tWe: {We}")
                    pass
        
        # Talley does not area-average local u_t; instead computes <u_t> with area-averaged parameters
        u_t_avg = 1.4 * np.cbrt(eps) * np.cbrt(6 * alpha_avg / ai_avg)
       # u_t_avg = self.area_avg('u_t_loc', method=avg_method)   #Quan, 1106
        We_avg = rho_f * u_t_avg**2 * (6 * alpha_avg / ai_avg) / sigma

        if u_t_avg > 0:

            sqrt_term = 1 - (We_cr / We_avg)  #Quan 1107
            if sqrt_term >= 0:
                COV_TI_avg = (u_t_avg * ai_avg**2 / alpha_avg) * np.sqrt(sqrt_term) * np.exp(-We_cr / We_avg)
            else:
                COV_TI_avg = 0  # or handle it in another way as needed

           # COV_TI_avg = (u_t_avg * ai_avg**2 / alpha_avg) * np.sqrt(1 - (We_cr / We_avg)) * np.exp(-We_cr / We_avg)

            COV_TI_loc_avg = self.area_avg('COV_TI_loc', method=avg_method)  #Quan, 1107
            if COV_TI_avg != 0 and not np.isnan(COV_TI_avg):
                COV_TI = COV_TI_loc_avg / COV_TI_avg
            else:
                COV_TI = 0  # or handle this case as needed

          #  COV_TI = self.area_avg('COV_TI_loc',method=avg_method) / COV_TI_avg
            
            if debug:
                print(f"\n\t\tu_t_avg: {u_t_avg}\tWe_avg: {We_avg}")
                print(f"\n\tCOV_TI: {COV_TI}\tCOV_TI_num: {self.area_avg('COV_TI_loc',method=avg_method)}\tCOV_TI_den: {COV_TI_avg}")
        else:
            COV_TI = 0

        self.COV_TI = COV_TI

        return COV_TI
    
    
    def reconstruct_void(self, method='talley', avg_method = 'legacy'):
        """ Reconstruct the void profile based on various methods

        Saves:
         * 'alpha_reconstructed' in midas_dict
         * self.roverRend
        
        Methods:
         * talley
         * double_linear (not_talley). Something I accidently invented while trying to implement Talley's method
        
        """

        debug = False

        if method.lower() == 'talley':
            self.roverRend = round(-1.472e-5 * self.Ref + 2.571,1)      # Inner r/R, outer end fixed at r/R = 1. Also, Talley rounds his r/R_end to the nearest 0.1

            # if debug:
            #     print(self.roverRend)

            def lineq(x, m, x0, b):
                return float(max(m * (x - x0) + b, 0))

            def find_alpha_peak(alpha_peak, rstar_peak=0.90):
                # Talley's reconstruction has a finer grid than experimental data
                # This will affect the slope of the drop-off point if r/R_end is not coincident with a point on the experiment mesh
                r_points = np.arange(0,1,0.05)
                self.add_mesh_points(r_points)

                interps = {}
                
                for angle, r_dict in self.data.items():                 # 360 degrees covered, not just one quadrant
                    
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
                        rstar_anchor = self.roverRend

                    elif angle_q1 == 0:
                        # For 0 degrees, value at r/R = 0 along the 90 degree axis
                        peak = lineq(x = 0, 
                                     m = 0 - alpha_peak / (self.roverRend - rstar_peak),
                                     x0 = rstar_peak,
                                     b = alpha_peak)

                        anchor = peak
                        rstar_anchor = 0

                    else:
                        # Find r* of previous angle at equivalent y-coordinate
                        y_peak = rstar_peak * np.sin(angle_q1 * np.pi / 180)
                        rstar_nn = y_peak / np.sin(angle_nn * np.pi / 180)

                        # if debug:
                        #     print(f"angle: {angle}\tangle_nn: {angle_nn}")

                        # Peak void fraction of current angle defined as void fraction at previous angle r*
                        peak = lineq(x = rstar_nn,
                                     m = m_nn,
                                     x0 = x0_nn,
                                     b = b_nn)
                        
                        if debug and angle == 22.5:
                            print(f"y_peak: {y_peak}\trstar_nn: {rstar_nn}")
                            pass
                        
                        anchor = 0
                        rstar_anchor = self.roverRend / np.cos((90 - angle_q1) * np.pi / 180)       # Talley's implementation in Excel
                        # rstar_anchor = self.roverRend / np.sin(angle_q1 * np.pi / 180)            # But why not like this
                        
                        # if self.roverRend > 0:
                        #     anchor = 0
                        #     rstar_anchor = self.roverRend
                        # else:
                        #     # Value at r/R = 0 along the 90 degree axis
                        #     anchor = lineq(x = 0,
                        #                    m = 0 - alpha_peak / (self.roverRend - rstar_peak),
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

                        if debug:
                            print(f"{angle}\t{rstar}:\t\talpha_rec: {midas_dict['alpha_reconstructed']:.4f}\talpha_dat: {midas_dict['alpha']:.4f}")
                    
                    self.alpha_peak = alpha_peak

                return abs( round(self.area_avg('alpha',method=avg_method),3) - self.area_avg('alpha_reconstructed',method=avg_method) )    # Talley's value is to three decimals of precision
            
            result = minimize(find_alpha_peak, x0 = 0.5, bounds = ((0,1),))      # The way they want me to format bounds is stupid. Python is stupid.

            if result.success:
                self.alpha_peak_reconstructed = result.x
                find_alpha_peak(self.alpha_peak_reconstructed)
            else:
                warnings.warn("Minimization did not return a successful result")
                print(result.message)
            
            print(f"\n\tr/R_end: {self.roverRend}")
            print(f"\talpha_peak: {self.alpha_peak}")
            print(f"\t⟨α⟩_data: {round(self.area_avg('alpha',method=avg_method),3)}")
            print(f"\t⟨α⟩_reconstructed: {self.area_avg('alpha_reconstructed',method=avg_method)}")

            '''
            def find_alpha_peak(alpha_peak):
                for angle, r_dict in self.data.items():
                    for rstar, midas_dict in r_dict.items():

                        x = rstar * np.cos(angle * np.pi / 180)
                        y = rstar * np.sin(angle * np.pi / 180)

                        # First calculate centerline void fraction (first linear interpolation)
                        if y > 0.9:
                            alpha_CL = alpha_peak / 0.1 * (1 - y)
                        else:
                            alpha_CL = max(alpha_peak / (0.9 - self.roverRend) * (y - self.roverRend), 0) # make sure it's not < 0. This covers for y < roverRend
                        alpha_CL = float(alpha_CL)

                        # Then calculate
                        if rstar >= 0.9:
                            xtrans = 0.9 * np.cos(angle * np.pi / 180)
                            midas_dict['alpha_reconstructed'] = alpha_CL / (1 - xtrans) * (1-x)
                        else:
                            midas_dict['alpha_reconstructed'] = alpha_CL

                return abs( self.area_avg('alpha') - self.area_avg('alpha_reconstructed') )
            
            result = minimize(find_alpha_peak, x0 = 0.5)

            if result.success:
                self.alpha_peak_reconstructed = result.x
                find_alpha_peak(self.alpha_peak_reconstructed)
            else:
                if debug:
                    warnings.warn("Minimization did not return a successful result")
                    print(f"⟨α⟩_data: {self.area_avg('alpha')}\n⟨α⟩_reconstructed: {self.area_avg('alpha_reconstructed')}\n")
            '''

        elif method.lower() == 'not_talley' or method.lower() == 'double_linear':
            self.roverRend = -1.472e-5 * self.Ref + 2.571

            def find_alpha_peak(alpha_peak):
                for angle, r_dict in self.data.items():
                    for rstar, midas_dict in r_dict.items():

                        x = rstar * np.cos(angle * np.pi / 180)
                        y = rstar * np.sin(angle * np.pi / 180)

                        # First calculate centerline void fraction (first linear interpolation)
                        if y > 0.9:
                            alpha_CL = alpha_peak / 0.1 * (1 - y)
                        else:
                            alpha_CL = max(alpha_peak / (0.9 - self.roverRend) * (y - self.roverRend), 0) # make sure it's not < 0. This covers for y < roverRend

                        # 
                        if np.sqrt(1 - y**2) == 0:
                            midas_dict['alpha_reconstructed'] = 0
                        else:
                            midas_dict['alpha_reconstructed'] = float(alpha_CL / np.sqrt(1 - y**2) * (np.sqrt(1 - y**2) - np.abs(x) ))

                return abs( self.area_avg('alpha') - self.area_avg('alpha_reconstructed') )
            
            result = minimize(find_alpha_peak, x0 = 0.5)

            if result.success:
                self.alpha_peak_reconstructed = result.x
                find_alpha_peak(self.alpha_peak_reconstructed)
            else:
                if debug:
                    warnings.warn("Minimization did not return a successful result")
                    print(f"⟨α⟩_data: {self.area_avg('alpha')}\n⟨α⟩_reconstructed: {self.area_avg('alpha_reconstructed')}\n")
            
        elif method.lower() == 'ryan':
            # TODO
            pass
        elif method.lower() == 'adix':
            
            def complicated_alpha(params):
                s, n, sigma = params
                for angle, r_dict in self.data.items():
                    for rstar, midas_dict in r_dict.items():

                        x = rstar * np.cos(angle * np.pi / 180)
                        y = rstar * np.sin(angle * np.pi / 180)
                        h = 1 - y
                        
                        if (np.sqrt(1-y**2) - abs(x)) < 0:
                            if debug: print((np.sqrt(1-y**2) - abs(x)), x, y)
                        # alpha = float(s * abs(np.sqrt(1-y**2) - abs(x))**(n) * np.sqrt(h) /sigma**2 * np.exp(-h**2 / (2*sigma**2)) )
                        alpha = float(abs(np.sqrt(1-y**2) - abs(x))**(n) *  s * 1 / sigma * (h/sigma)**(0.5) * np.exp(-(h/sigma)**0.5))
                        if alpha < 1e-3:
                            alpha = 0
                        elif alpha > 1:
                            alpha = 1
                        midas_dict['alpha_reconstructed'] = alpha

                self.calc_errors('alpha', 'alpha_reconstructed')

                # return (abs( self.area_avg('alpha') - self.area_avg('alpha_reconstructed') )/self.area_avg('alpha') + abs(self.max('alpha') - self.max('alpha_reconstructed')) / self.max('alpha')) * 100
                return self.sum('eps_sq_alpha_alpha_reconstructed')

                
            # result = minimize(complicated_alpha, x0 = [0.2, 2, 0.3], method='Nelder-Mead', bounds=( (0.01, 1), (1, 100), (0.01, 10) ), options = {'maxiter': 100000, 'disp': True}, tol = 1e-9)
            result = minimize(complicated_alpha, x0 = [0.041, 0.1417, 0.0155], method='Nelder-Mead', bounds=( (0.0001, 1), (0.001, 10), (0.0001, 0.1) ), options = {'maxiter': 100000})

            if result.success:
                self.reconstruct_s = result.x[0]
                self.reconstruct_n = result.x[1]
                self.reconstruct_sigma = result.x[2]
                complicated_alpha(result.x)
            else:
                if debug:
                    warnings.warn("Minimization did not return a successful result")
                    complicated_alpha(result.x)
                    print(f"⟨α⟩_data: {self.area_avg('alpha')}\n⟨α⟩_reconstructed: {self.area_avg('alpha_reconstructed')}\n{result.x}")


        return self.area_avg("alpha_reconstructed")
    
    def plot_profiles2(self, param, save_dir = '.', show=True, x_axis='vals', errorbars = 0.0, 
                      const_to_plot = [90, 67.5, 45, 22.5, 0], include_complement = True, 
                      fig_size=(4,4), title=True, label_str = '', legend_loc = 'best', xlabel_loc = 'center', include_const = False,
                      set_min = None, set_max = None, show_spines = True, xlabel_loc_coords = None, ylabel_loc_coords = None, cs=None) -> None:
        """ Simplified plot_profiles with no rotation option

        Inputs:
         - param, any parameter string in midas_dict, e.g. 'alpha', 'ai', etc., or a list of parameter strings ['alpha_G1', 'alpha_G2']
        
        Options:
         - save_dir, directory in which to save the .png file. Will not save the file unless show = False
         - show, display the figure (in an iPython notebook or have it pop up)
         - x_axis, the variable to put on the x-axis. Usually for vertical flow this is 'rs', for horizontal 'vals'. Also can be 'phis'
         - const_to_plot, a list of angles (if x-axis is 'vals' or 'rs') or rs (if x-axis is 'phis')
         - include_complement, includes the complementary angle of the const_to_plot (i.e. 270 with 90)
         - fig_size, a tuple of the figure size, in inches
         - title, puts a title above the figure
         - label_str, a string that will replace the default title
         - legend_loc, passed to ax.legend
         - xlable_loc, passed to ax.xaxis.set_label_coords
         - xlable_loc_coords, passed to ax.xaxis.set_label_coords
         - errorbars, value of errorbars to put on the data
         
        
        """
        # TODO rewrite so it always loops over a list of params. If there's only one, just put it in a list at the begininng 

        plt.rcParams.update({'font.size': 10})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"

        fig, ax = plt.subplots(figsize=fig_size, dpi=300, layout='compressed')

        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        # Tick marks facing in
        ax.tick_params(direction='in',which='both')

        if errorbars > 0:
            ax.plot([], [], ' ', label = f"{errorbars*100:0.1f}% error bars") # dummy to just get this text in the legend

        ms = marker_cycle()
        if cs is None or type(param) == list:
            cs = color_cycle()
        elif cs == 'infer' and type(param) == str:
            cs = color_cycle(set_color = param)

        if set_min == None and type(param) == str:
            set_min = self.min(param)
        elif set_min == None and type(param) == list:
            mins = []
            for specific_param in param:
                mins.append(self.min(specific_param))

            set_min = min(mins)
        
        if set_max == None and type(param) == str:
            set_max = self.max(param) *1.1

        elif set_min == None and type(param) == list:
            maxs = []
            for specific_param in param:
                maxs.append(self.max(specific_param))

            set_max = max(maxs) *1.1

        if x_axis == 'vals' or x_axis == 'r':
            for angle in const_to_plot:
                r_dict = self.data[angle]
                rs = []
                vals = []
                for r, midas_output in r_dict.items():
                    rs.append(r)

                    if type(param) == str:
                        try:
                            vals.append(midas_output[param])
                        except:
                            if abs(r - 1) < 0.0001:
                                vals.append(0.0)
                            else:
                                vals.append(0.0)
                                print(f"Could not find {param} for φ = {angle}, r = {r}. Substituting 0")
                    elif type(param) == list:
                        for i, specific_param in enumerate(param):
                            vals.append([])
                            try:
                                vals[i].append(midas_output[specific_param])
                            except:
                                if abs(r - 1) < 0.0001:
                                    vals[i].append(0.0)
                                else:
                                    vals[i].append(0.0)
                                    print(f"Could not find {specific_param} for φ = {angle}, r = {r}. Substituting 0")

                if include_complement:
                    if angle > 180:
                        print("Error: Cannot find complement to angle > 180. Skipping")
                    else:
                        r_dict = self.data[angle+180]
                        for r, midas_output in r_dict.items():
                            rs.append(-r)
                            if type(param) == str:
                                try:
                                    vals.append(midas_output[param])
                                except:
                                    if abs(r - 1) < 0.0001:
                                        vals.append(0.0)
                                    else:
                                        vals.append(0.0)
                                        print(f"Could not find {param} for φ = {angle}, r = {r}. Substituting 0")
                            elif type(param) == list:
                                for i, specific_param in enumerate(param):
                                    vals.append([])
                                    try:
                                        vals[i].append(midas_output[specific_param])
                                    except:
                                        if abs(r - 1) < 0.0001:
                                            vals[i].append(0.0)
                                        else:
                                            vals[i].append(0.0)
                                            print(f"Could not find {specific_param} for φ = {angle}, r = {r}. Substituting 0")
                if type(param) == str:
                    vals = [var for _, var in sorted(zip(rs, vals))]
                    rs = sorted(rs)
                    errs = errorbars * np.abs(np.asarray(vals))

                    if x_axis == 'vals':
                        ax.errorbar(vals, rs, xerr = errs, capsize=3, ecolor = "black", label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
                    elif x_axis == 'r':
                        ax.errorbar(rs, vals, yerr = errs, capsize=3, ecolor = "black", label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
                
                elif type(param) == list:
                    temp = []
                    for vals_list in vals:
                        temp_list = [var for _, var in sorted(zip(rs, vals_list))]
                        temp.append(temp_list)
                    rs = sorted(rs)

                    for specific_param, val_list in zip(param, temp):
                        errs = errorbars * np.abs(np.asarray(val_list))
                        # print(val_list)
                        if cs == 'infer':
                            cs = color_cycle(set_color = specific_param)
                        
                        
                        legend_str = ''
                        if 'alpha' in specific_param:
                            specific_param = specific_param.replace('alpha', r'$\alpha$')
                        elif 'ai' in specific_param:
                            specific_param = specific_param.replace('ai', r'$a_{i}$')
                        elif 'ug' in specific_param:
                            specific_param = specific_param.replace('ug', r'$v_{g}$')
                        elif 'vf' in specific_param:
                            specific_param = specific_param.replace('vf', r'$v_{f}$')
                        
                        elif specific_param == 'vr':
                            legend_str = r'$v_{r, G1}$'
                        elif specific_param == 'vr2':
                            legend_str = r'$v_{r, G2}$'
                        elif specific_param == 'Dsm1':
                            legend_str = r'$D_{sm1}$'
                        elif specific_param == 'Dsm2':
                            legend_str = r'$D_{sm2}$'

                        if legend_str == '':
                            legend_str = specific_param

                        if include_const:
                            legend_str += f", {angle}°"

                        if x_axis == 'vals':
                            ax.errorbar(val_list, rs, xerr = errs, capsize=3, ecolor = "black", color=next(cs), marker=next(ms), linestyle = '--', label = legend_str)
                        elif x_axis == 'r':
                            ax.errorbar(rs, val_list, yerr = errs, capsize=3, ecolor = "black", color=next(cs), marker=next(ms), linestyle = '--', label = legend_str)
                    
            
            if x_axis == 'vals':
                ax.set_ylim(-1, 1)
                ax.set_xlim(set_min, set_max)
            elif x_axis == 'r':
                ax.set_xlim(-1, 1)
                ax.set_ylim(set_min, set_max)
            
        
        elif x_axis == 'phi':
            ax.plot([], [], label=r'$r/R$', color='white', linestyle = None)
            for rtarget in const_to_plot:
                phis = []
                vals = []
                for angle, r_dict in self.data.items():
                    for rstar, midas_output in r_dict.items():
                        if abs(rstar - rtarget) < 0.001:
                            phis.append(angle)
                            try:
                                vals.append(midas_output[param])
                            except:
                                if abs(rstar - 1) < 0.0001:
                                    vals.append(0.0)
                                else:
                                    vals.append(0.0)
                                    print(f"Could not find {param} for φ = {angle}, r = {rstar}. Substituting 0")

                vals = [var for _, var in sorted(zip(phis, vals))]
                phis = sorted(phis)
                
                ax.plot(phis, vals, label=f'{rtarget:0.1f}', color=next(cs), marker=next(ms), linestyle = '--')
            
            ax.set_xlim(0, 360)
            ax.set_ylim(set_min, set_max)
                
        else:
            print(f"invalid axis for plot_profiles: {x_axis}. Current supported options are 'rs', 'vals' and 'phi'")
            return
        
        if label_str == '' and type(param) == str:
            if param == 'alpha':
                label_str = r'$\alpha\ [-]$'
            elif param == 'ai':
                label_str = r'$a_{i}\ [1/m]$'
            elif param == 'ug1' or param == 'ug2':
                label_str = r'$v_{g}\ [m/s]$'
            elif param == 'vf':
                label_str = r'$v_{f}\ [m/s]$'
            elif param == 'vr' or param == 'vr2':
                label_str = r'$v_{r}\ [m/s]$'
            elif param == 'Dsm1':
                label_str = r'$D_{sm1}\ [mm]$'
            else:
                label_str = param
        
        if x_axis == 'vals':
            ax.set_xlabel(label_str, loc = xlabel_loc)
            ax.set_ylabel(r'$r/R$ [-]')
            ax.set_yticks(np.arange(-1, 1.01, 0.2))
            
            ax.spines['bottom'].set_position(('data', 0))
            
            if set_min == 0 or set_max == 0:
                ax.spines['left'].set_position(('data', set_min))
                ax.spines['right'].set_position(('data', 0))
                # ax.yaxis.tick_right()
                # ax.yaxis.set_label_position("right")

            if not show_spines:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
            else:
                ax2 = ax.twiny()
                ax2.get_xaxis().set_visible(False)


        elif x_axis == 'r':
            ax.set_ylabel(label_str, loc = xlabel_loc)
            ax.set_xlabel(r'$r/R$ [-]')
            ax.set_xticks(np.arange(-1, 1.01, 0.2))

            ax.spines['bottom'].set_position(('data', max(0, set_min)))
            
            ax.spines['left'].set_position(('data', 0))
            
            if not show_spines:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
            else:
                ax2 = ax.twiny()
                ax2.get_xaxis().set_visible(False)
            
            
        elif x_axis == 'phi':
            
            ax.set_ylabel(label_str, loc = xlabel_loc)
            ax.set_xlabel(r'$\varphi$ [-]')

            ax.set_xticks([0, 90, 180, 270, 360])
               
        if hasattr(title,'__len__'):
            # Set your own title! -DHK
            ax.set_title(title)
        elif title:
            ax.set_title(self.name)

        
        

        leg = ax.legend(loc=legend_loc, edgecolor='white')

        if xlabel_loc_coords:
            ax.xaxis.set_label_coords(*xlabel_loc_coords)

        if ylabel_loc_coords:
            ax.yaxis.set_label_coords(*ylabel_loc_coords)


        
        ax.set_aspect('auto', adjustable='datalim', share=True)
        #fake_ax.set_box_aspect(1)
        
        if show:
            plt.show()
        else:
            if type(param) == str:
                plt.savefig(os.path.join(save_dir, f'{param}_profile_vs_{x_axis}_{self.name}.png'))
            else:
                plt.savefig(os.path.join(save_dir, f'{"_".join(param)}_profile_vs_{x_axis}_{self.name}.png'))
            plt.close()
        return

    def plot_profiles(self, param, save_dir = '.', show=True, x_axis='vals', 
                      const_to_plot = [90, 67.5, 45, 22.5, 0], include_complement = True, 
                      rotate=False, fig_size=(4,4), title=True, label_str = '', legend_loc = 'best', xlabel_loc = 'center',
                      set_min = None, set_max = None, show_spines = True, force_RH_y_axis = False, xlabel_loc_coords = None, cs=None) -> None:
        """ Plot profiles of param over x_axis, for const_to_plot, i.e. α over r/R for φ = [90, 67.5 ... 0]. 
        
        Include_complement will continue with the negative side if x_axis = 'r' 

        Also has an option to rotate the graph based on self.theta, but it's a little sketchy
        
        """

        self.mirror()
        plt.rcParams.update({'font.size': 10})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"

        log_x = False # This breaks, so I removed it from the arguments to the function

        

        if rotate:

            # Set up the figure to be rotated by theta
            import matplotlib as mpl
            from matplotlib.transforms import Affine2D
            import mpl_toolkits.axisartist.floating_axes as floating_axes
            fig = plt.figure(figsize=fig_size)
            if x_axis == 'r':
                plot_extents = self.min(param), self.max(param)*1.1, -1, 1
                transform = Affine2D().scale(fig_size[0] / (self.max(param)*1.1 - self.min(param)), fig_size[1] / (1 - -1)).rotate_deg(self.theta)
            else:
                plot_extents = self.min(param), self.max(param)*1.1, 0, 360
                transform = Affine2D().scale(fig_size[0] / (self.max(param)*1.1 - self.min(param)), fig_size[1] / (360-0)).rotate_deg(self.theta)
            
            
            helper = floating_axes.GridHelperCurveLinear(transform, plot_extents)
            fake_ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=helper)
            ax = fake_ax.get_aux_axes(transform)

        else:
            fig, ax = plt.subplots(figsize=fig_size, dpi=300, layout='compressed')
            fake_ax = ax

        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        # Tick marks facing in
        ax.tick_params(direction='in',which='both')

        ms = marker_cycle()
        if cs is None:
            cs = color_cycle()
        elif cs == 'infer':
            cs = color_cycle(set_color = param)
        else:
            print("I hope cs is a generator that returns valid colors")

        if set_min == None:
            set_min = self.min(param)
        
        if set_max == None:
            set_max = self.max(param) *1.1

        if x_axis == 'vals' or x_axis == 'r':
            for angle in const_to_plot:
                r_dict = self.data[angle]
                rs = []
                vals = []
                for r, midas_output in r_dict.items():
                    rs.append(r)
                    try:
                        vals.append(midas_output[param])
                    except:
                        if abs(r - 1) < 0.0001:
                            vals.append(0.0)
                        else:
                            vals.append(0.0)
                            print(f"Could not find {param} for φ = {angle}, r = {r}. Substituting 0")

                if include_complement:
                    if angle > 180:
                        print("Error: Cannot find complement to angle > 180. Skipping")
                    else:
                        r_dict = self.data[angle+180]
                        for r, midas_output in r_dict.items():
                            rs.append(-r)
                            try:
                                vals.append(midas_output[param])
                            except:
                                if abs(r - 1) < 0.0001:
                                    vals.append(0.0)
                                else:
                                    vals.append(0.0)
                                    print(f"Could not find {param} for φ = {angle}, r = {r}. Substituting 0")

                vals = [var for _, var in sorted(zip(rs, vals))]
                rs = sorted(rs)
                if x_axis == 'vals':
                    ax.plot(vals, rs, label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
                if x_axis == 'r':
                    ax.plot(rs, vals, label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
            
            if x_axis == 'vals':
                ax.set_ylim(-1, 1)
                ax.set_xlim(set_min, set_max)
            elif x_axis == 'r':
                ax.set_xlim(-1, 1)
                ax.set_ylim(set_min, set_max)
                fake_ax.set_xlim(-1, 1)
                fake_ax.set_ylim(set_min, set_max)
            
        
        elif x_axis == 'phi':
            ax.plot([], [], label=r'$r/R$', color='white', linestyle = None)
            for rtarget in const_to_plot:
                phis = []
                vals = []
                for angle, r_dict in self.data.items():
                    for rstar, midas_output in r_dict.items():
                        if abs(rstar - rtarget) < 0.001:
                            phis.append(angle)
                            try:
                                vals.append(midas_output[param])
                            except:
                                if abs(rstar - 1) < 0.0001:
                                    vals.append(0.0)
                                else:
                                    vals.append(0.0)
                                    print(f"Could not find {param} for φ = {angle}, r = {rstar}. Substituting 0")

                vals = [var for _, var in sorted(zip(phis, vals))]
                phis = sorted(phis)
                if rotate:
                    ax.plot(vals, phis, label=f'{rtarget:0.1f}', color=next(cs), marker=next(ms), linestyle = '--')
                else:
                    ax.plot(phis, vals, label=f'{rtarget:0.1f}', color=next(cs), marker=next(ms), linestyle = '--')
            if rotate:
                ax.set_ylim(0, 360)
                ax.set_xlim(set_min, set_max)
                
            else:
                ax.set_xlim(0, 360)
                ax.set_ylim(set_min, set_max)
                
        else:
            print(f"invalid axis for plot_profiles: {x_axis}. Current supported options are 'r' and 'phi'")
            return
        
        if label_str == '':
            if param == 'alpha':
                label_str = r'$\alpha\ [-]$'
            elif param == 'ai':
                label_str = r'$a_{i}\ [1/m]$'
            elif param == 'ug1':
                label_str = r'$v_{g}\ [m/s]$'
            else:
                label_str = param
        
        if x_axis == 'vals':
            fake_ax.set_xlabel(label_str, loc = xlabel_loc)
            fake_ax.set_ylabel(r'$r/R$ [-]')
            fake_ax.set_yticks(np.arange(-1, 1.01, 0.2))
            #fake_ax.set_xticks(np.linspace(self.min(param), self.max(param), 7))

        elif x_axis == 'r':
            fake_ax.set_ylabel(label_str, loc = xlabel_loc)
            fake_ax.set_xlabel(r'$r/R$ [-]')
            fake_ax.set_xticks(np.arange(-1, 1.01, 0.2))

        elif x_axis == 'phi':
            if not rotate:
                fake_ax.set_ylabel(label_str, loc = xlabel_loc)
                fake_ax.set_xlabel(r'$\varphi$ [-]')

                fake_ax.set_xticks([0, 90, 180, 270, 360])
                #fake_ax.set_yticks(np.linspace(self.min(param), self.max(param), 7))
            else:
                fake_ax.set_xlabel(label_str, loc = xlabel_loc)
                fake_ax.set_ylabel(r'$\varphi$ [-]')
                
                fake_ax.set_yticks([0, 90, 180, 270, 360])
                # fake_ax.set_xticks(np.linspace(self.min(param), self.max(param), 7))

                #fake_ax.tick_params(axis='both', labelrotation=-self.theta)
        
        if title:
            ax.set_title(self.name)

        ax.spines['bottom'].set_position(('data', 0))


        if set_min == 0 or set_max == 0 or force_RH_y_axis:
            ax.spines['left'].set_position(('data', set_min))
            ax.spines['right'].set_position(('data', 0))
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position("right")

        if not show_spines:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        else:
            ax2 = ax.twiny()
            ax2.get_xaxis().set_visible(False)


        if log_x:
            ax.set_xlim(self.min(param, nonzero=True), self.max(param)*1.2)
            ax.set_xscale('log')
            fake_ax.set_xscale('log')
        
        if rotate:
            fig.add_subplot(fake_ax)
        ax.legend(loc=legend_loc, edgecolor='white')

        if xlabel_loc_coords:
            ax.xaxis.set_label_coords(*xlabel_loc_coords)

        fake_ax.set_aspect('auto', adjustable='datalim', share=True)
        ax.set_aspect('auto', adjustable='datalim', share=True)
        #fake_ax.set_box_aspect(1)
        
        if show:
            plt.show()
        else:
            plt.savefig(os.path.join(save_dir, f'{param}_profile_vs_{x_axis}_{self.name}.png'))
            plt.close()
        return   

    def plot_isoline(self, param:str, iso_axis:str, iso_val:float, fig_size=4, plot_res=100, 
                     save_dir = '.', show=True, extra_text = '') -> None:
        """ Plot profiles of param over iso_axis at iso_val 
        
        Based on interpolation, so plot_res changes the resolution
        
        """

        fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300, layout='compressed')

        plt.rcParams.update({'font.size': 10})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"

        self.mirror()

        if iso_axis == 'x':
            lim = np.sqrt(1 - iso_val**2)
            ys = np.linspace(-lim, lim, plot_res)
            xs = iso_val * np.ones(ys.size)
            plt.plot(ys, self(xs, ys, param, interp_method = 'linear_xy'), color = self.marker_color, marker=None, linestyle= self.line_style)

            plt.xlim(-1, 1)
            plt.xlabel(r'$y/R$ [-]')
            plt.ylabel(param)
            plt.show()

        elif iso_axis == 'y':
            lim = np.sqrt(1 - iso_val**2)
            xs = np.linspace(-lim, lim, plot_res)
            ys = iso_val * np.ones(xs.size)
            plt.plot(xs, self(xs, ys, param, interp_method = 'linear_xy'), color = self.marker_color, marker=None, linestyle= self.line_style)

            plt.xlim(-1, 1)
            plt.xlabel(r'$x/R$ [-]')
            plt.ylabel(param)
            plt.show()

        if show:
            plt.show()
        else:
            plt.savefig( os.path.join(save_dir, f"{param}_iso_{iso_axis}={iso_val}_{self.name + extra_text}.png") )
            plt.close()
        return

    def plot_contour(self, param:str, save_dir = '.', show=True, set_max = None, set_min = None, fig_size = 4, label_str = None, suppress_colorbar = False,
                     rot_angle = 0, ngridr = 50, ngridphi = 50, colormap = 'hot_r', num_levels = 0, level_step = 0.01, title = False, title_str = '', extra_text = '',
                     annotate_h = False, cartesian = False, h_star_kwargs = {'method': 'max_dsm', 'min_void': '0.05'}, plot_measured_points = False, font_size = 16) -> None:
        
        """Method to plot contour of a given param
        
        Generates a contour plot of any parameter in midas_dict, e.g. 'alpha', 'ai', etc. By default, just shows the figure,
        but if a save_dir is specified, it will save it there instead. label_str can adjust the label of the colorbar. Can accept Latex format, e.g.
        r"$\alpha$ [-]"

        set_max and set_min set the bounds of the contour plot, and the colormap option allows for any colors that matplotlib supports. ngridr, ngridphi,
        num_levels, all adjust how fine the contour plot is generated.

        annotate_h is an option to draw a horizontal line at some given position. This was implemented when investigating where the bubble layer typically
        stops. h_star_kwargs and cartesian are optinos related to this.

        plot_measured points is neat, it plots circles where original data was detected (determined prior to mirroring)
        
        """ 
        if cartesian:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300)
        else:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300, subplot_kw=dict(projection='polar'))
        plt.rcParams.update({'font.size': font_size})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"
        
        # self.mirror()
        rs = []
        phis = []
        vals = []

        for phi_angle, r_dict in self.data.items():
            for r, midas_output in r_dict.items():
                if r >= 0:
                    rs.append(r)
                    phis.append(phi_angle)
                    try:
                        vals.append(midas_output[param])
                    except:
                        if abs(r - 1) < 0.0001:
                            vals.append(0.0)
                        else:
                            vals.append(0.0)
                            print(f"Could not find {param} for φ = {phi_angle}, r = {r}. Substituting 0")

        rs = np.asarray(rs)
        phis = (np.asarray(phis) + rot_angle) * np.pi / 180 
        vals = np.asarray(vals)   

        ri = np.linspace(0, 1, ngridr)
        phii = np.linspace(rot_angle * np.pi / 180, 2*np.pi + rot_angle * np.pi / 180, ngridphi)

        triang = tri.Triangulation(phis, rs)
        interpolator = tri.LinearTriInterpolator( triang, vals )

        PHII, RI = np.meshgrid(phii, ri)
        XI = RI * np.cos(PHII)
        YI = RI * np.sin(PHII)
        parami = interpolator(PHII, RI)

        if debug and False: 
            print(self.name, f'{param} to contour plot', file=debugFID)
            print(Xi, Yi, alphai, file=debugFID)
            print(f'max {param}', np.amax(alphai), file = debugFID)
        
        if param == 'alpha' and (np.nanmax(vals) > 1.0 or np.nanmax(vals) > 1.0):
            np.savetxt(f'contour_error_dump_{self.name}.csv', vals, delimiter=',')
            #print(vals)
            raise(ValueError('Alpha exceeds 1.0!\nSaving problematic array'))
            print('Alpha exceeds 1.0!\nSaving problematic array')
            return
        
        extend_min = False
        extend_max = False
        
        if set_min == None:
            set_min = np.min(parami)

        if set_max == None:
            set_max = np.max(parami) + (np.max(parami) * 0.1)

        if set_min > np.min(parami):
            extend_min = True

        if set_max < np.max(parami):
            extend_max = True

        if extend_max and extend_min:
            extend_opt = 'both'
        elif extend_min and not extend_max and set_min > np.min(parami):
            extend_opt = 'min'
        elif extend_max and not extend_min:
            extend_opt = 'max'
        elif not extend_min and not extend_max:
            extend_opt = 'neither'

        if abs(set_max - set_min) < level_step:
            print(f'Warning: Level step {level_step} too larger for range {set_max}, {set_min}. Defaulting to 0.01*(set_max - set_min)')
            level_step = 0.01*(set_max - set_min)

        if num_levels:
            lvs = np.linspace(set_min, set_max, num_levels)
        else:
            if (abs(set_max) < 1e-8): # if the max is 0, start the counting from there
                lvs = np.arange(set_max, abs(set_min) + 1e-8, level_step)
                lvs = np.flip(lvs) * -1
            else:
                lvs = np.arange(set_min, set_max + 1e-8, level_step)
        
        if cartesian:
            mpbl = ax.contourf(XI, YI, parami, levels = lvs, vmin = set_min, vmax = set_max, cmap = colormap)

            x = np.linspace(-1, 1, 100)
            # Make a circle to look nice
            plt.plot(x, np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
            plt.plot(x, -np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
        else:
            mpbl = ax.contourf(PHII, RI, parami, levels = lvs, 
                               vmin = set_min, vmax = set_max, cmap = colormap, extend=extend_opt)
            if plot_measured_points: 
                # print(np.asarray(self.original_mesh)[:,0]* np.pi/180 , np.asarray(self.original_mesh)[:,1])
                measured_thetas = np.asarray(self.original_mesh)[:,0]* np.pi/180
                measured_rs = np.asarray(self.original_mesh)[:,1]

                for i, r in enumerate(np.asarray(self.original_mesh)[:,1]):
                    if r < 0:
                        measured_rs[i] = - measured_rs[i]
                        measured_thetas[i] = (measured_thetas[i] + np.pi) % (2*np.pi)

                ax.plot( measured_thetas, measured_rs, marker='o', fillstyle = 'none', mec = 'red', linewidth = 0, ms = 2)

        if annotate_h:
            if not cartesian:
                print("Warning: annotate_h assumes that we're plotting in Cartesian. Annotation may not be correct")
            hstar = self.find_hstar_pos(**h_star_kwargs)
            ax.annotate('h', 
                    (-1,1-hstar), 
                    (1,1-hstar),
                    arrowprops=dict(color='r', width=3, headwidth=3),
                    color='r',
                    verticalalignment='center'
        )
            
        ax.grid(False)
        if cartesian:
            plt.axis('square')
            plt.xlabel (r'$x/R$ [-]')
            plt.ylabel(r'$y/R$ [-]')
            ax.grid(False)
        else:
            ax.set_yticklabels([])
            ax.set_xticklabels([])

        #plt.clim(vmin=set_min, vmax=set_max)
        
        if label_str == None:
            label_str = param

        if not suppress_colorbar:
            tx_step = round((set_max - set_min)/5,-int(np.floor(np.log10((set_max - set_min)/10))))
            tx = np.arange(set_min,set_max,tx_step)

            fig.colorbar(mpbl, label=label_str, ticks=tx)

        if title_str != '':
            title = True
        
        if title:
            if title_str == '':
                plt.title(self.name)
            else:
                plt.title(title_str)

        plt.tight_layout()
        
        #cb = plt.colorbar(ticks = [0, 0.05, 0.1, 0.15, 0.2])

        if show:
            plt.show()
        else:
            plt.savefig( os.path.join(save_dir, f"{param}_contours_{self.name + extra_text}.png") )
            plt.close()
        return ax

    def plot_surface(self, param:str, save_dir = '.', show=True, set_max = None, set_min = None, rotate_gif=False, elev_angle = 145, 
                     azim_angle = 0, roll_angle = 180, title=True, ngridr = 50, ngridphi = 50, 
                     plot_surface_kwargs = None, solid_color = False, label_str = None, title_str = '') -> None:
        """Method to plot a surface of a given param
        
        Can save a static image or rotating gif, starting at elev_angle, azim_angle, roll_angle. These angles also 
        the viewing angle for the static image.

        Can specify a label or title str.

        Plot_surface_kwargs is how to specify vmin, vmax, colormap, etc.
        
        """ 
        
        if plot_surface_kwargs is None:
            plot_surface_kwargs = {}
        plt.rcParams.update({'font.size': 16})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"

        self.mirror()
        rs = []
        phis = []
        vals = []

        for phi_angle, r_dict in self.data.items():
            for r, midas_output in r_dict.items():
                if r >= 0:
                    rs.append(r)
                    phis.append(phi_angle)

                    try:
                        vals.append(midas_output[param])
                    except:
                        vals.append(np.NaN)
                        print(f"Could not find {param} for φ = {phi_angle}, r = {r}. Substituting NaN")

        rs = np.asarray(rs)
        phis = (np.asarray(phis)) * np.pi / 180 
        vals = np.asarray(vals)   

        ri = np.linspace(0, 1, ngridr)
        phii = np.linspace(0, 2*np.pi, ngridphi)

        triang = tri.Triangulation(phis, rs)
        interpolator = tri.LinearTriInterpolator( triang, vals )

        PHII, RI = np.meshgrid(phii, ri)
        parami = interpolator(PHII, RI)

        Xi = RI * np.cos(PHII)
        Yi = RI * np.sin(PHII)

        fig, ax = plt.subplots(figsize=(5, 5), subplot_kw={"projection": "3d"})

        if 'vmin' not in plot_surface_kwargs.keys(): 
            plot_surface_kwargs.update({'vmin': self.min(param)})

        if 'vmax' not in plot_surface_kwargs.keys(): 
            plot_surface_kwargs.update({'vmax': self.max(param)})

        if 'cmap' not in plot_surface_kwargs.keys() and not solid_color: 
            plot_surface_kwargs.update({'cmap': 'viridis'})

        surf = ax.plot_surface(Xi, Yi, parami, **plot_surface_kwargs)
        
        #plt.legend()
        ax.set_xlabel (r'$x/R$ [-]')
        ax.set_ylabel(r'$y/R$ [-]')

        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        
        ax.set_zlim([plot_surface_kwargs['vmin'], plot_surface_kwargs['vmax']])
        
        if set_min == None:
            set_min = np.min(parami)

        if set_max == None:
            set_max = np.max(parami) + (np.max(parami) * 0.1)

        tx_step = round((set_max - set_min)/5,-int(np.floor(np.log10((set_max - set_min)/10))))
        tx = np.arange(set_min,set_max,tx_step)
        
        if label_str:
            ax.set_zlabel(label_str)
            fig.colorbar(surf, label=label_str,ticks=tx)
        else:
            ax.set_zlabel(param)
            fig.colorbar(surf, label=param,ticks=tx)
        
        if title: 
            if title_str:
                plt.title(title_str)
            else:
                plt.title(self.name)

        ax.view_init(azim=azim_angle, roll=roll_angle, elev=elev_angle)
               
        #cb = plt.colorbar(ticks = [0, 0.05, 0.1, 0.15, 0.2])
        if show:
            plt.show()
        else:
            plt.savefig( os.path.join(save_dir, f"{param}_surface_{self.name}.png") )
        
        if rotate_gif:
            import matplotlib.animation as animation
            
            def rotate(angle):
                ax.view_init(azim=angle)
            
            rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
            rot_animation.save(os.path.join(save_dir, f'{param}_surface_rotation_{self.name}.gif'), dpi=80)

        return

    def plot_spline_contour(self, param:str, save_dir = '.', show=True, set_max = None, set_min = None, fig_size = 4,
                     rot_angle = 0, ngridr = 50, ngridphi = 50, colormap = 'hot_r', num_levels = 100, title = False,
                     annotate_h = False, cartesian = False, h_star_kwargs = {'method': 'max_dsm', 'min_void': '0.05'},
                     grad = 'None') -> None:
        """Plots a contour from a spline interpolation

        Will fit the spline if necessary. 
        
        """ 
        
        if param not in self.spline_interp.keys():
            print(f"Warning: {param} not found in spline_interp dict, running fit_spline")
            self.fit_spline(param)

        phii = np.linspace(0, 2*np.pi, 100)
        phii_arg = phii + rot_angle * np.pi/180
        ri = np.linspace(0, 1, 100)

        if cartesian:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300)
        else:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300, subplot_kw=dict(projection='polar'))
        plt.rcParams.update({'font.size': 10})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"

        PHII, RI = np.meshgrid(phii, ri)
        XI = RI * np.cos(PHII)
        YI = RI * np.sin(PHII)

        if grad == 'None':
            VALS = (self.spline_interp[param](phii_arg, ri)).T
        elif grad == 'r':
            VALS = (self.spline_interp[param](phii_arg, ri, dy=1)).T
        elif grad == 'r':
            VALS = (self.spline_interp[param](phii_arg, ri, dy=1)).T
        elif grad == 'phi':
            VALS =  (1./RI *self.spline_interp[param](phii_arg, ri, dx = 1)).T
        elif grad == 'phinor':
            VALS =  (self.spline_interp[param](phii_arg, ri, dx = 1)).T
        elif grad == 'y':
            VALS = (self.spline_interp['alpha'](phii_arg, ri, dx=1) * np.cos(phii_arg)*ri / (ri**2+1e-8) + self.spline_interp['alpha'](phii_arg, ri, dy=1) * np.sin(phii_arg)).T
        else:
            print(f"Error: unrecognized grad type {grad}")


        if cartesian:
            plt.contourf(XI, YI, VALS, levels = num_levels, vmin = set_min, vmax = set_max, cmap = colormap)

            x = np.linspace(-1, 1, 100)
            # Make a circle to look nice
            plt.plot(x, np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
            plt.plot(x, -np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
        else:
            plt.contourf(PHII, RI, VALS, levels = num_levels, vmin = set_min, vmax = set_max, cmap = colormap)

        if annotate_h:
            if not cartesian:
                print("Warning: annotate_h assumes that we're plotting in Cartesian. Annotation may not be correct")
            hstar = self.find_hstar_pos(**h_star_kwargs)
            ax.annotate('h', 
                    (-1,1-hstar), 
                    (1,1-hstar),
                    arrowprops=dict(color='r', width=3, headwidth=3),
                    color='r',
                    verticalalignment='center'
        )
            
        ax.grid(False)
        if cartesian:
            plt.axis('square')
            plt.xlabel (r'$x/R$ [-]')
            plt.ylabel(r'$y/R$ [-]')
            ax.grid(True)
        else:
            ax.set_yticklabels([])
            ax.set_xticklabels([])

        if grad == 'None':
            param_label = param
        else:
            param_label = "grad_" + param + "_" + grad 
        plt.colorbar(label=param_label)
        
        if title:
            plt.title(self.name)

        plt.tight_layout()
        
        #cb = plt.colorbar(ticks = [0, 0.05, 0.1, 0.15, 0.2])

        if show:
            plt.show()
        else:
            plt.savefig( os.path.join(save_dir, f"{param_label}_spline_contours_{self.name}.png") )
            plt.close()
        return

    def rough_FR_ID(self) -> None:
        """Identifies the flow regime for the given condition, by some rough methods

        Inputs:
         - None

        Stores:
         - "FR" in self

        Returns:
         - FR

        First checks if it matches any given by previous researchers, or the hierarchical
        clustering algorithm results

        1 = bubbly
        2 = plug
        3 = slug
        4 = churn
        5 = stratified
        6 = stratified wavy
        7 = annular
        
        """ 

        # These values are taken from the respective theses where the data comes from
               
        drew_frid = [
            (0, 4,    0.111, 1),
            (0, 4,    0.186, 1),
            (0, 4,    0.333, 1),
            (0, 1.56, 0.111, 2),
            (0, 2.5,  0.111, 2),
            (0, 3.5,  0.33,  3),
            (0, 2,    0.33,  2),
            (0, 2,    1.0,   3),
            (0, 1,    0.33,  2),
            
            (30, 4,    0.111, 1),
            (30, 4,    0.186, 1),
            (30, 4,    0.333, 1),
            (30, 1.56, 0.111, 2),
            (30, 2.5,  0.111, 2),
            (30, 3.5,  0.33,  1),
            (30, 2,    0.33,  3),
            (30, 2,    1.0,   3),
            (30, 1,    0.33,  3),

            (60, 4,    0.111, 1),
            (60, 4,    0.186, 1),
            (60, 4,    0.333, 1),
            (60, 1.56, 0.111, 3),
            (60, 2.5,  0.111, 2), #I changed this, because the point was a consistent outlier
            (60, 3.5,  0.33,  1),
            (60, 2,    0.33,  3),
            (60, 2,    1.0,   3),
            (60, 1,    0.33,  3),

            (80, 4,    0.111, 1),
            (80, 4,    0.186, 1),
            (80, 4,    0.333, 1),
            (80, 1.56, 0.111, 1),
            (80, 2.5,  0.111, 1),
            (80, 3.5,  0.33,  1),
            (80, 2,    0.33,  3),
            (80, 2,    1.0,   3),
            (80, 1,    0.33,  2),

            (90, 4,    0.111, 1),
            (90, 4,    0.186, 1),
            (90, 4,    0.333, 1),
            (90, 1.56, 0.111, 1),
            (90, 2.5,  0.111, 1),
            (90, 3.5,  0.33,  1),
            (90, 2,    0.33,  1),
            (90, 2,    1.0,   1),
            (90, 1,    0.33,  1)
        ]

        kong_frid = [
            (0, 4.0, 0.111, 1),
            (0, 4.0, 0.25,  2),
            (0, 5.0, 0.327, 1),
            (0, 6.0, 1.5,   1),
            (0, 6.0, 0.297, 1),
            (0, 4.0, 0.575, 2),
            (0, 3.0, 0.189, 2),
            (0, 3.0, 0.332, 2),
            (0, 2.5, 0.329, 2),
            (0, 4.0, 1.040, 3),
            (0, 3.5, 0.799, 3),
            (0, 3.5, 1.030, 3),

            (0, 4.00, 0.111, 1),
            (0, 4.00, 0.186, 1),
            (0, 5.00, 0.105, 1),
            (0, 5.00, 0.167, 1),
            (0, 5.00, 0.327, 1),
            (0, 6.00, 0.089, 1),
            (0, 6.00, 0.149, 1),
            (0, 6.00, 0.297, 1),
            (0, 6.00, 0.575, 1),
            (0, 5.00, 0.575, 1),
            (0, 6.00, 1.000, 1),
            (0, 5.00, 1.000, 1),
            (0, 4.00, 0.330, 1),
            (0, 3.50, 0.332, 1),
            (0, 3.00, 0.327, 2),
            (0, 4.00, 0.507, 1),
            (0, 4.00, 1.000, 3),
            (0, 4.00, 3.000, 3),
            (0, 2.00, 0.330, 2),
            (0, 2.00, 0.575, 2),
            (0, 2.00, 1.000, 3),
            (0, 2.00, 3.000, 3),
            (0, 1.00, 0.330, 2)
        ]

        talley_frid = [
            (0, 6, 0.15, 1), 
            (0, 6, 0.25, 1), 
            (0, 6, 0.5,  1), 
            (0, 6, 1.0,  1), 
            (0, 5, 0.15, 1), 
            (0, 5, 0.25, 1), 
            (0, 5, 0.50, 1), 
            (0, 4, 0.15, 1), 
            (0, 4, 0.25, 1)
        ]

        # This one is taken from David's hierarchical clustering algorithm results
        frid = [
            (0 , 0.04, 0.02  ,5),  
            (0 , 0.04, 0.07  ,5),  
            (0 , 0.04, 0.24  ,5),  
            (0 , 0.04, 0.76  ,5),  
            (0 , 0.04, 2.44  ,6),  
            (0 , 0.04, 7.81  ,6),  
            (0 , 0.04, 20.7  ,7),  
            (0 , 0.05, 7.81  ,6),  
            (0 , 0.05, 20.4  ,7),  
            (0 , 0.07, 7.81  ,6),  
            (0 , 0.07, 19.6  ,7),  
            (0 , 0.08, 0.02  ,5),  
            (0 , 0.08, 0.07  ,5),  
            (0 , 0.08, 0.24  ,2),  
            (0 , 0.08, 0.76  ,5),  
            (0 , 0.08, 2.44  ,6),  
            (0 , 0.08, 7.81  ,6),  
            (0 , 0.08, 19.3  ,7),  
            (0 , 0.10, 7.81  ,6),  
            (0 , 0.10, 18.6  ,7),  
            (0 , 0.13, 7.81  ,6),  
            (0 , 0.13, 18.2  ,7),  
            (0 , 0.17, 0.02  ,5),  
            (0 , 0.17, 0.07  ,2),  
            (0 , 0.17, 0.24  ,2),  
            (0 , 0.17, 0.76  ,2),  
            (0 , 0.17, 2.44  ,6),  
            (0 , 0.17, 17.6  ,7),  
            (0 , 0.18, 7.81  ,6),  
            (0 , 0.38, 0.02  ,2),  
            (0 , 0.38, 0.07  ,2),  
            (0 , 0.38, 0.24  ,2),  
            (0 , 0.38, 0.76  ,3),  
            (0 , 0.38, 2.44  ,3),  
            (0 , 0.38, 7.81  ,3),  
            (0 , 0.38, 14.7  ,3),  
            (0 , 0.83, 0.02  ,2),  
            (0 , 0.83, 0.07  ,2),  
            (0 , 0.83, 0.24  ,2),  
            (0 , 0.83, 0.76  ,3),  
            (0 , 0.83, 2.44  ,3),  
            (0 , 1.82, 0.02  ,1),  
            (0 , 1.82, 0.07  ,2),  
            (0 , 1.82, 0.24  ,2),  
            (0 , 1.82, 0.76  ,2),  
            (0 , 1.82, 2.44  ,3),  
            (0 , 2.50, 0.11  ,3),  
            (0 , 4.00, 0.02  ,1),  
            (0 , 4.00, 0.07  ,1),  
            (0 , 4.00, 0.111 ,1),  
            (0 , 4.00, 0.186 ,1),  
            (0 , 4.00, 0.24  ,1),  
            (0 , 4.00, 0.33  ,1),  
            (30, 0.04, 0.02  ,2),  
            (30, 0.04, 0.07  ,3),  
            (30, 0.04, 0.24  ,3),  
            (30, 0.04, 0.76  ,3),  
            (30, 0.04, 7.81  ,3),  
            (30, 0.04, 21.2  ,7),  
            (30, 0.05, 7.81  ,3),  
            (30, 0.05, 22.5  ,7),  
            (30, 0.07, 7.81  ,3),  
            (30, 0.07, 21.7  ,7),  
            (30, 0.08, 0.02  ,2),  
            (30, 0.08, 0.07  ,2),  
            (30, 0.08, 0.24  ,3),  
            (30, 0.08, 0.76  ,3),  
            (30, 0.08, 2.44  ,3),  
            (30, 0.08, 7.81  ,3),  
            (30, 0.08, 19.9  ,7),  
            (30, 0.10, 7.81  ,3),  
            (30, 0.10, 20.9  ,7),  
            (30, 0.13, 7.81  ,3),  
            (30, 0.13, 20.2  ,7),  
            (30, 0.17, 0.02  ,2),  
            (30, 0.17, 0.07  ,2),  
            (30, 0.17, 0.24  ,3),  
            (30, 0.17, 0.76  ,3),  
            (30, 0.17, 2.44  ,3),  
            (30, 0.17, 7.81  ,3),  
            (30, 0.17, 25.0  ,7),  
            (30, 0.38, 0.07  ,2),  
            (30, 0.38, 0.24  ,2),  
            (30, 0.38, 0.76  ,3),  
            (30, 0.38, 2.44  ,3),  
            (30, 0.38, 7.81  ,3),  
            (30, 0.38, 15.5  ,3),  
            (30, 0.83, 0.02  ,1),  
            (30, 0.83, 0.07  ,2),  
            (30, 0.83, 0.24  ,2),  
            (30, 0.83, 0.76  ,3),  
            (30, 0.83, 2.44  ,3),  
            (30, 1.82, 0.02  ,1),  
            (30, 1.82, 0.07  ,1),  
            (30, 1.82, 0.24  ,2),  
            (30, 1.82, 0.76  ,3),  
            (30, 1.82, 2.44  ,3),  
            (30, 4.00, 0.02  ,1),  
            (30, 4.00, 0.07  ,1),  
            (30, 4.00, 0.24  ,1),  				  
            (60, 0.04, 0.02  ,3),  
            (60, 0.04, 0.07  ,3),  
            (60, 0.04, 0.24  ,3),  
            (60, 0.04, 0.76  ,3),  
            (60, 0.04, 2.44  ,3),  
            (60, 0.04, 7.81  ,3),  
            (60, 0.04, 21.2  ,7),  
            (60, 0.05, 7.81  ,3),  
            (60, 0.05, 20.9  ,7),  
            (60, 0.07, 6.72  ,3),  
            (60, 0.07, 20.1  ,7),  
            (60, 0.08, 0.02  ,3),  
            (60, 0.08, 0.07  ,3),  
            (60, 0.08, 0.24  ,3),  
            (60, 0.08, 0.76  ,3),  
            (60, 0.08, 2.44  ,3),  
            (60, 0.08, 6.73  ,3),  
            (60, 0.08, 19.3  ,7),  
            (60, 0.10, 6.75  ,3),  
            (60, 0.10, 19.8  ,7),  
            (60, 0.13, 6.76  ,3),  
            (60, 0.13, 18.2  ,7),  
            (60, 0.17, 0.02  ,3),  
            (60, 0.17, 0.07  ,3),  
            (60, 0.17, 0.24  ,3),  
            (60, 0.17, 0.76  ,3),  
            (60, 0.17, 2.12  ,3),  
            (60, 0.17, 6.78  ,3),  
            (60, 0.17, 16.7  ,7),  
            (60, 0.38, 0.02  ,3),  
            (60, 0.38, 0.07  ,3),  
            (60, 0.38, 0.24  ,3),  
            (60, 0.38, 0.76  ,3),  
            (60, 0.38, 2.14  ,3),  
            (60, 0.38, 7.81  ,3),  
            (60, 0.38, 13.5  ,7),  
            (60, 0.83, 0.02  ,1),  
            (60, 0.83, 0.07  ,1),  
            (60, 0.83, 0.24  ,3),  
            (60, 0.83, 0.76  ,3),  
            (60, 0.83, 2.44  ,3),  
            (60, 1.82, 0.02  ,1),  
            (60, 1.82, 0.07  ,1),  
            (60, 1.82, 0.24  ,1),  
            (60, 1.82, 0.76  ,3),  
            (60, 1.82, 2.44  ,3),  
            (60, 4.00, 0.02  ,1),  
            (60, 4.00, 0.07  ,1),  
            (60, 4.00, 0.24  ,1),  
            (80, 0.04, 0.02  ,1),  
            (80, 0.04, 0.07  ,3),  
            (80, 0.04, 0.24  ,3),  
            (80, 0.04, 0.76  ,3),  
            (80, 0.04, 2.44  ,4),  
            (80, 0.04, 7.81  ,4),  
            (80, 0.04, 26.4  ,7),  
            (80, 0.05, 7.81  ,4),  
            (80, 0.05, 23.8  ,7),  
            (80, 0.07, 7.81  ,4),  
            (80, 0.07, 23.1  ,7),  
            (80, 0.08, 0.02  ,1),  
            (80, 0.08, 0.07  ,3),  
            (80, 0.08, 0.24  ,3),  
            (80, 0.08, 0.76  ,3),  
            (80, 0.08, 2.44  ,4),  
            (80, 0.08, 7.81  ,4),  
            (80, 0.08, 22.4  ,7),  
            (80, 0.10, 7.81  ,4),  
            (80, 0.10, 20.6  ,7),  
            (80, 0.13, 7.81  ,4),  
            (80, 0.13, 20.1  ,7),  
            (80, 0.17, 0.02  ,1),  
            (80, 0.17, 0.07  ,3),  
            (80, 0.17, 0.24  ,3),  
            (80, 0.17, 0.76  ,3),  
            (80, 0.17, 2.44  ,4),  
            (80, 0.17, 7.81  ,4),  
            (80, 0.17, 19.1  ,7),  
            (80, 0.38, 0.02  ,1),  
            (80, 0.38, 0.07  ,1),  
            (80, 0.38, 0.24  ,1),  
            (80, 0.38, 0.76  ,3),  
            (80, 0.38, 2.44  ,3),  
            (80, 0.38, 7.81  ,3),  
            (80, 0.38, 15.1  ,4),  
            (80, 0.83, 0.02  ,1),  
            (80, 0.83, 0.07  ,1),  
            (80, 0.83, 0.24  ,1),  
            (80, 0.83, 0.76  ,3),  
            (80, 0.83, 2.44  ,3),  
            (80, 1.82, 0.02  ,1),  
            (80, 1.82, 0.07  ,1),  
            (80, 1.82, 0.24  ,1),  
            (80, 1.82, 0.76  ,3),  
            (80, 1.82, 2.44  ,3),  
            (80, 4.00, 0.02  ,1),  
            (80, 4.00, 0.07  ,1),  
            (80, 4.00, 0.24  ,1),  
            (90, 0.04, 0.02  ,1),  
            (90, 0.04, 0.07  ,3),  
            (90, 0.04, 0.24  ,3),  
            (90, 0.04, 0.76  ,3),  
            (90, 0.04, 2.44  ,4),  
            (90, 0.04, 7.81  ,4),  
            (90, 0.04, 21.5  ,7),  
            (90, 0.08, 0.02  ,1),  
            (90, 0.08, 0.07  ,3),  
            (90, 0.08, 0.24  ,3),  
            (90, 0.08, 0.76  ,3),  
            (90, 0.08, 2.44  ,4),  
            (90, 0.08, 7.81  ,4),  
            (90, 0.08, 20.2  ,7),  
            (90, 0.17, 0.02  ,1),  
            (90, 0.17, 0.07  ,3),  
            (90, 0.17, 0.24  ,3),  
            (90, 0.17, 0.76  ,3),  
            (90, 0.17, 2.44  ,4),  
            (90, 0.17, 7.81  ,4),  
            (90, 0.17, 18.2  ,7),  
            (90, 0.38, 0.02  ,1),  
            (90, 0.38, 0.07  ,1),  
            (90, 0.38, 0.24  ,1),  
            (90, 0.38, 0.76  ,3),  
            (90, 0.38, 2.44  ,3),  
            (90, 0.38, 7.81  ,3),  
            (90, 0.83, 0.02  ,1),  
            (90, 0.83, 0.07  ,1),  
            (90, 0.83, 0.24  ,1),  
            (90, 0.83, 0.76  ,3),  
            (90, 0.83, 2.44  ,3),  
            (90, 1.82, 0.02  ,1),  
            (90, 1.82, 0.07  ,1),  
            (90, 1.82, 0.24  ,1),  
            (90, 1.82, 0.76  ,3),  
            (90, 4.00, 0.02  ,1),  
            (90, 4.00, 0.07  ,1),  
            (90, 4.00, 0.24  ,1),  
            (0 , 0.04, 0.04  ,5),  
            (0 , 0.04, 0.13  ,5),
            (0 , 0.04, 0.43  ,5),
            (0 , 0.04, 1.36  ,5),
            (0 , 0.04, 4.37  ,6),  
            (0 , 0.04, 13.98 ,6),
            (0 , 0.05, 0.04  ,5),  
            (0 , 0.05, 0.13  ,5),  
            (0 , 0.05, 0.43  ,2),  
            (0 , 0.05, 1.36  ,3),  
            (0 , 0.05, 4.37  ,6),
            (0 , 0.05, 13.98 ,6),  
            (0 , 0.07, 0.04  ,5),   
            (0 , 0.07, 0.13  ,5),  
            (0 , 0.07, 0.43  ,2),  
            (0 , 0.07, 1.36  ,3),  
            (0 , 0.07, 4.37  ,6),  
            (0 , 0.07, 13.98 ,7),  
            (0 , 0.08, 0.04  ,5),  
            (0 , 0.08, 0.13  ,5),  
            (0 , 0.08, 0.43  ,2),  
            (0 , 0.08, 1.36  ,3),  
            (0 , 0.08, 4.37  ,6),  
            (0 , 0.08, 13.98 ,7),  
            (0 , 0.17, 0.04  ,2),  
            (0 , 0.38, 0.04  ,2),  
            (0 , 0.83, 0.04  ,2),  
            (0 , 1.82, 0.04  ,1),  
            (0 , 4.00, 0.04  ,1),  
            (0 , 4.00, 0.111 ,1),  
            (0 , 4.00, 0.186 ,1),  
            (0 , 4.00, 0.33  ,1),  
            (0 , 2.50, 0.111 ,1),  
            (0 , 1.56, 0.111 ,2),  
            (0 , 0.08, 0.76  ,2),  
            (0 , 0.08, 0.07  ,5),  
            (0 , 0.04, 0.24  ,2),  
            (30, 0.04, 0.04  ,2),  
            (30, 0.04, 0.13  ,3),  
            (30, 0.04, 0.43  ,3),  
            (30, 0.04, 1.36  ,3),  
            (30, 0.04, 4.37  ,3),  
            (30, 0.04, 13.98 ,7),  
            (30, 0.05, 0.04  ,2),  
            (30, 0.05, 0.13  ,3),  
            (30, 0.05, 0.43  ,3),  
            (30, 0.05, 1.36  ,3),  
            (30, 0.05, 4.37  ,3),  
            (30, 0.05, 13.98 ,7),  
            (30, 0.07, 0.04  ,2),  
            (30, 0.07, 0.13  ,3),  
            (30, 0.07, 0.43  ,3),  
            (30, 0.07, 1.36  ,3),  
            (30, 0.07, 4.37  ,3),  
            (30, 0.07, 13.98 ,7),  
            (30, 0.08, 0.04  ,2),  
            (30, 0.08, 0.13  ,3),  
            (30, 0.08, 0.43  ,3),  
            (30, 0.08, 1.36  ,3),  
            (30, 0.08, 4.37  ,3),  
            (30, 0.08, 13.98 ,7),  
            (30, 0.17, 0.04  ,2),  
            (30, 0.38, 0.04  ,2),  
            (30, 0.83, 0.04  ,2),  
            (30, 1.82, 0.04  ,1),  
            (30, 4.00, 0.04  ,1),  
            (30, 4.00, 0.111 ,1),  
            (30, 4.00, 0.186 ,1),  
            (30, 4.00, 0.33  ,1),  
            (30, 2.50, 0.111 ,1),  
            (30, 1.56, 0.111 ,2),  
            (30, 0.38, 0.02  ,2),  
            (30, 0.38, 0.76  ,3),  
            (30, 0.17, 7.81  ,3),  
            (30, 0.17, 17.99 ,7),  
            (30, 0.04, 2.44  ,3),  
            (60, 0.04, 0.04  ,3),  
            (60, 0.04, 0.13  ,3),  
            (60, 0.04, 0.43  ,3),  
            (60, 0.04, 1.36  ,3),  
            (60, 0.04, 4.37  ,3),  
            (60, 0.04, 13.98 ,7),  
            (60, 0.05, 0.04  ,3),  
            (60, 0.05, 0.13  ,3),  
            (60, 0.05, 0.43  ,3),  
            (60, 0.05, 1.36  ,3),  
            (60, 0.05, 4.37  ,3),  
            (60, 0.05, 13.98 ,7),  
            (60, 0.07, 0.04  ,3),  
            (60, 0.07, 0.13  ,3),  
            (60, 0.07, 0.43  ,3),  
            (60, 0.07, 1.36  ,3),  
            (60, 0.07, 4.37  ,3),  
            (60, 0.07, 13.98 ,7),  
            (60, 0.08, 0.04  ,3),  
            (60, 0.08, 0.13  ,3),  
            (60, 0.08, 0.43  ,3),  
            (60, 0.08, 1.36  ,3),  
            (60, 0.08, 4.37  ,3),  
            (60, 0.08, 13.98 ,7),  
            (60, 0.17, 0.04  ,3),  
            (60, 0.38, 0.04  ,3),  
            (60, 0.83, 0.04  ,1),  
            (60, 1.82, 0.04  ,1),  
            (60, 4.00, 0.04  ,1),  
            (60, 4.00, 0.111 ,1),  
            (60, 4.00, 0.186 ,1),  
            (60, 4.00, 0.33  ,1),  
            (60, 2.50, 0.111 ,1),  
            (60, 1.56, 0.111 ,1),  
            (80, 0.04, 0.04  ,3),  
            (80, 0.04, 0.13  ,3),  
            (80, 0.04, 0.43  ,3),  
            (80, 0.04, 1.36  ,4),  
            (80, 0.04, 4.37  ,4),  
            (80, 0.04, 13.98 ,7),  
            (80, 0.05, 0.04  ,1),  
            (80, 0.05, 0.13  ,3),  
            (80, 0.05, 0.43  ,3),  
            (80, 0.05, 1.36  ,3),  
            (80, 0.05, 4.37  ,4),  
            (80, 0.05, 13.98 ,7),  
            (80, 0.07, 0.04  ,1),  
            (80, 0.07, 0.13  ,3),  
            (80, 0.07, 0.43  ,3),  
            (80, 0.07, 1.36  ,3),  
            (80, 0.07, 4.37  ,4),  
            (80, 0.07, 13.98 ,7),  
            (80, 0.08, 0.04  ,1),  
            (80, 0.08, 0.13  ,3),  
            (80, 0.08, 0.43  ,3),  
            (80, 0.08, 1.36  ,3),  
            (80, 0.08, 4.37  ,4),  
            (80, 0.08, 13.98 ,7),  
            (80, 0.17, 0.04  ,1),  
            (80, 0.38, 0.04  ,1),  
            (80, 0.83, 0.04  ,1),  
            (80, 1.82, 0.04  ,1),  
            (80, 4.00, 0.04  ,1),  
            (80, 4.00, 0.111 ,1),  
            (80, 4.00, 0.186 ,1),  
            (80, 4.00, 0.33  ,1),  
            (80, 2.50, 0.111 ,1),  
            (80, 1.56, 0.111 ,1),
            (90, 0.04, 0.048 ,1),  
            (90, 0.04, 0.158 ,3),  
            (90, 0.04, 0.531 ,3),  
            (90, 0.04, 1.719 ,3),  
            (90, 0.04, 5.541 ,4),  
            (90, 0.04, 17.66 ,7),  
            (90, 0.05, 0.047 ,1),  
            (90, 0.05, 0.159 ,3),  
            (90, 0.05, 0.538 ,3),  
            (90, 0.05, 1.729 ,3),  
            (90, 0.05, 5.577 ,4),  
            (90, 0.05, 17.57 ,7),  
            (90, 0.07, 0.047 ,1),  
            (90, 0.07, 0.161 ,3),  
            (90, 0.07, 0.537 ,3),  
            (90, 0.07, 1.724 ,3),  
            (90, 0.07, 5.562 ,4),  
            (90, 0.07, 17.39 ,7),  
            (90, 0.08, 0.047 ,1),  
            (90, 0.08, 0.158 ,3),  
            (90, 0.08, 0.536 ,3),  
            (90, 0.08, 1.725 ,3),  
            (90, 0.08, 5.564 ,4),  
            (90, 0.08, 17.36 ,7),  
            (90, 0.17, 0.047 ,1),  
            (90, 0.38, 0.051 ,1),  
            (90, 0.83, 0.051 ,1),  
            (90, 1.82, 0.048 ,1),  
            (90, 4.00, 0.046 ,1),  
            (90, 4.00, 0.127 ,1),  
            (90, 4.00, 0.209 ,1),  
            (90, 4.00, 0.372 ,1),  
            (90, 2.50, 0.131 ,1),  
            (90, 1.56, 0.136 ,1),  
            (90, 0.13, 20.74 ,7),  
            (90, 0.13, 9.870 ,4),  
            (90, 0.10, 20.50 ,7),  
            (90, 0.10, 9.899 ,4),  
            (90, 0.07, 21.75 ,7),  
            (90, 0.07, 9.916 ,4),  
            (90, 0.05, 22.60 ,7),  
            (90, 0.05, 9.941 ,4),  
            (90, 4.00, 0.111 ,1),  
            (90, 4.00, 0.186 ,1),  
            (90, 4.00, 0.33  ,1),  
            (90, 2.50, 0.111 ,1),  
            (90, 1.56, 0.111 ,1) 
        ]


        for theta_ID, jf_ID, jg_ID, FR_ID in drew_frid:
            if self.theta == theta_ID:
                if abs( self.jf - jf_ID) < 0.05 and abs(self.jgloc - jg_ID) < 0.05:
                    self.FR = FR_ID
                    return

        for theta_ID, jf_ID, jg_ID, FR_ID in kong_frid:
            if self.theta == theta_ID:
                if abs( self.jf - jf_ID) < 0.05 and abs(self.jgloc - jg_ID) < 0.05:
                    self.FR = FR_ID
                    return

        for theta_ID, jf_ID, jg_ID, FR_ID in talley_frid:
            if self.theta == theta_ID:
                if abs( self.jf - jf_ID) < 0.05 and abs(self.jgloc - jg_ID) < 0.05:
                    self.FR = FR_ID
                    return
                
        for theta_ID, jf_ID, jg_ID, FR_ID in frid:
            if self.theta == theta_ID:
                if abs( self.jf - jf_ID) < 0.05 and abs(self.jgloc - jg_ID) < 0.05:
                    self.FR = FR_ID
                    return


        if self.theta == 90:
            self.FR = 1
        else:
            self.FR = 2

        if (self.jf > 3.4) and (self.jgloc < 1.0):
            self.FR = 1
            if self.theta == 0 and self.jf < 3.7:
                self.FR = 2
        
        #print(f"Warning: Flow regime not accurately identified for {self.name}, defaulting to {self.FR}")
        return self.FR

    def TD_FR_ID(self) -> None:
        """ This function estimates the flow regime based on Taitel and Dukler's method

        Not implemented yet
        
        """
        #dpdxL

        return self.FR




def print_tab_keys() -> None:
    """Convenience function in case you forget a tab key
    """

    print(tab_keys)
    return

tab_keys = [
    'roverR',
    'time',
    'frequency',
    'num_spherical',
    'num_distorted',
    'num_cap',
    'num_slug',
    'num_G1',
    'num_G2',
    'num_total',
    'obs_0',
    'obs_1',
    'obs_2',
    'obs_3',
    'bub_freq',
    'pair_spherical',
    'pair_distorted',
    'pair_cap',
    'pair_slug',
    'total_paired',
    'percent_paired',
    'alpha_spherical',
    'alpha_distorted',
    'alpha_cap',
    'alpha_slug',
    'alpha_G1',
    'alpha_G2',
    'alpha',
    'ai_spherical',
    'ai_distorted',
    'ai_cap',
    'ai_slug',
    'ai_G1',
    'ai_G2',
    'ai',
    'Dsm1',
    'Lcl1',
    'Dsm2',
    'Lcl2',
    'ug1',
    'ug2',
    'sigma_ug1',
    'sigma_ug2',
    'fluctuation',
    'alpha_ug1',
    'alpha_ug2',
    'alpha_ug',
    'alpha_Dsm1',
    'alpha_Dsm2',
    'r01',
    'r02',
    'r03',
    'r12',
    'r13',
    'r23',
    'vf',
    'jf_loc',
    'jf',
    'delta_p',
    'sigma_delta_p',
    'vr'
]

old_tab_keys = [
    'roverR',
    'time',
    'num_spherical',
    'num_distorted',
    'num_cap',
    'num_total',
    'bub_freq',
    'obs_1',
    'obs_2',
    'obs_3',
    'total_paired',
    'percent_paired',
    'alpha_spherical',
    'alpha_distorted',
    'alpha_cap',
    'alpha',
    'ai_spherical',
    'ai_distorted',
    'ai_cap',
    'ai',
    'Dsm1',
    'Dsm2',
    'ug1',
    'ug2',
    'sigma_ug1',
    'sigma_ug2',
    'fluctuation_ug',
]

zero_data = dict(zip(tab_keys, [0]*len(tab_keys)))
