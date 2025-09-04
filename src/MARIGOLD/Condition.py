from .config import *
from .operations import *

class Condition:
    """
    Class to handle the local probe data

    Data is stored in the ``Condition.data`` property. It's actually 3 layers of dictionary

    ``self.data[angle]`` gives a dictionary with the various r/R

    ``self.data[angle][r/R]`` gives a dictionary with the MIDAS output (``midas_dict``)

    The MIDAS output is itself a dictionary, with the keys listed in the "tab_keys" array
    So data[angle][r/R]['alpha'] should give you the void fraction at r/R for phi = angle
    This structure is initialized with zeros for the MIDAS output at the pipe center and wall

    Can also get the data at a local point from calling the condition, syntax

    ``self(phi, r, 'param')``. Phi is in radians, the arguments can be constants or numpy arrays.
    Also has an option for interpolation, 'interp_method' 
    """

    debugFID = None
    
    def __init__(self, jgref:float, jgloc:float, jf:float, theta:int, port:str, database:str, fluids = 'air-water', g = 9.81, p_atm = 101325, T = 293.15) -> None:
        """Initialize condition object
        
        **Args**:
        
         - ``jgref``: reference superficial gas velocity.
         - ``jgloc``: local superficial gas velocity.
         - ``jf``: superficial liquid velocity.
         - ``theta``: facility angle
         - ``port``: string to represent the port the data was taken at.
         - ``database``: database that the data was a part of
         - ``fluids``: fluids to use. Defaults to 'air-water'.
         - ``g``: gravitational acceleration. Defaults to 9.81.
        
        **Raises**:
        
         - ``NotImplementedError``: Invalid fluid pair
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

        if database == 'Ryan' or database == 'ryan' or database == 'adix' or database == 'neup':
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

        # Constants
        if fluids == 'air-water' or fluids == 'water-air':
            self.rho_f = 998        # kg/m^3
            self.rho_g = 1.204      # kg/m^3
            self.mu_f = 0.001002    # Pa-s
            self.mu_g = 0.01803e-3  # Pa-s
            self.R_spec = 287.058   # J/kg-K
            self.sigma = 0.0728     # N/m

            self.Ref = self.rho_f * self.jf * self.Dh / self.mu_f
            self.Reg = self.rho_g * self.jgloc * self.Dh / self.mu_g
        
        else:
            raise NotImplementedError(f"{fluids} not available, try 'air-water'")
        
        self.g = g
        self.gz = g * np.sin(np.radians(self.theta))
        self.p_atm = p_atm
        self.T = T

        return

    def __eq__(self, __o: object) -> bool:
        if isinstance(__o, Condition):

            return ((self.jf == __o.jf) and (self.jgref == __o.jgref) and (self.theta == __o.theta) and (self.port == __o.port) and self.database == __o.database)
        
        return False

    def __hash__(self) -> int:
        """hash
        
        
        **Returns**:
        
         - hash(repr(self))
        """
        return hash(repr(self))

    def __repr__(self) -> str:
        """representation
        
        
        **Returns**:
        
         - ``self.name``, "jf=X.XX_jgloc=Y.YY_theta=ZZ_port=PX_Database"
        """
        return self.name

    def __call__(self, phi_in:np.ndarray, r_in:np.ndarray, param:str, interp_method='None') -> np.ndarray:
        """Returns the value of param at :math:`(\\varphi, r)`. Phi is in radians, r nondimensional
        
        **Args**:
        
         - ``phi_in``: phi values to return values of
         - ``r_in``: r values to return values of
         - ``param``: ``midas_dict`` parameter to return values of. See :any:`print_params` for options
         - ``interp_method``: method to interpolate :math:`(\\varphi, r)`. Defaults to 'None'.

             - ``'None'``, just use the data
             - ``'linear'``, linear interpolation
             - ``'spline'``, spline interpolation
             - ``'linear_xy'``, phi = x, r = y
        
        **Raises**:
        
         - ``NameError``: if invalid ``interp_method`` selected.
        
        **Returns**:
        
         - Values of ``param`` 
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
        """_Prints out all the information in a Condition in a structured way 
        
        **Args**:
        
         - ``print_to_file``: option to print the output to a file. Defaults to False.
         - ``FID``: file ID to write to. Defaults to debugFID.
         - ``mirror``: mirror before printing. Defaults to False.
        """
        # Prints out all the information in a Condition in a structured way 

        # Specifically, everything in the Condition.phi dictionary, which has angles
        # and r/Rs.

        # Can either print to a file (specified by FID) or to stdout. Option to mirror the 
        # data, if that hasn't already been done

        # 

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
        """Mirrors data, so we have data for every angle.

        Typically, data is only recorded in one or two quadrants of the pipe cross section. This function copies that data across specified lines of symmetry to ensure data is present at every angle
        
        **Args**:
        
         - ``method``: method, or lines of symmetry to assume. Defaults to ``'sym90'``. Options:
             
             - ``'sym90'``, for horizontal or inclined flows which exhibit symmetry across the :math:`\\varphi = 90\\degree` line 
             - ``'axisym'``, for vertical upward or downwards axisymmetric flows
             - ``'avg_axisym'``, for vertical upward or downwards axisymmetric flows. Averages the angle with the most data and its complement.

         - ``sym90``: overwrites method, for backwards compatibility. Defaults to False.
         - ``axisym``: overwrites method, for backwards compatibility. Defaults to False.
         - ``uniform_rmesh``: ensure every angle has data for every r/R point. Will linearly interpolate when data on either side is available. This only  considers the +r/R mesh. Defaults to False.
         - ``uniform_rmesh_fill``: what to fill in for ``uniform_rmesh``. Defaults to 'interp'. Could also be 0.
         - ``force_remirror``: this overrides the typical behavior of exiting if ``mirror`` has been called before. Defaults to False.
        """
        #  Mirrors data, so we have data for every angle

        # :param method: method to use for mirroring. Options are 'sym90', 'axisym', and 'avg_axisym'. Defaults to 'sym90'
        # :type method: str, optional
        # :param sym90: overwrites method, for backwards compatibility, defaults to False
        # :type sym90: bool, optional
        # :param axisym: overwrites method, for backwards compatibility, defaults to False
        # :type axisym: bool, optional
        # :param uniform_rmesh: ensure every angle has data for every r/R point. Will linearly interpolate when data on either side is available. This only  considers the +r/R mesh, defaults to False
        # :type uniform_rmesh: bool, optional
        # :param force_remirror: Usually can only mirror once. This overrides that restriction, defaults to False
        # :type force_remirror: bool, optional

        # Saves:
        #  - original_mesh
        # 

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

        # Clean up None's
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                for param, value in midas_dict.items():
                    if value is None:
                        print(f"Warning: 'None' value for {param} at {rstar}. Setting to 0")
                        midas_dict.update( {param : 0})

        return
    
    def add_mesh_points(self, r_points:list, suppress=False):
        """Method for adding additional r/R points

        Data linearly interpolated based on surrounding data (specifically using :any:`__call__` at the :math:`(\\varphi, r)` location in question)
        
        **Args**:
        
         - ``r_points``: List of r locations to add
         - ``suppress``: print KeyErrors. Defaults to False.
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
    
    def fit_spline(self, param: str) -> None:
        """fits a spline to the data of param using ``scipy.interpolate.RectBivariateSpline``
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to fit spline. See :func:`~MARIGOLD.Condition.print_params` for options
        
        **Raises**:
        
         - ``KeyError``: If invalid parameter called
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
        """Fit a linear interpolant to the experimental data
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
         - ``suppress``: Suppress Key Error. Defaults to ``True``.
         - ``refit``: if ``True``, will fit the linear interpolant, even if it exists already. Defaults to ``False``.
        
        **Raises**:
        
         - ``KeyError``: If ``param`` not available
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
        """fits a linear Cartesian interpreter to the data
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
         - ``suppress``: Suppress ``KeyErrors``. Defaults to True.
        
        **Raises**:
        
         - ``KeyError``: if the requested ``param`` is not present everywhere
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
        """Sums values of ``param``
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to sum. See :func:`~MARIGOLD.Condition.print_params` for options
        
        **Returns**:
        
         - sum
        """
        total = 0
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                total += midas_dict[param]

        return total
    
    def max(self, param: str, recalc=True) -> float:
        """Finds the maximum value of ``param``
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to find maximum of. See :func:`~MARIGOLD.Condition.print_params` for options
         - ``recalc``: recalculate maximum. If False, and this has been called before, returns previously calculated value. Defaults to True.
        
        **Returns**:
        
         - Maximum value
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
        """Finds the location of the maximum value of ``param``
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to find maximum of. See :func:`~MARIGOLD.Condition.print_params` for options
        
        **Returns**:
        
         - ``(rstar, angle)`` location corresponding to the maximum
        """
        max = 0
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] > max:
                    max = midas_dict[param]
                    location = (rstar, angle)

        return (location)
    
    def min(self, param: str, recalc = True, nonzero=False)-> float:
        """Finds the minimum value of ``param``
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to find minimum of. See :func:`~MARIGOLD.Condition.print_params` for options
         - ``recalc``: recalculate minimum. If False, and this has been called before, returns previously calculated value. Defaults to True.
         - ``nonzero``: If True, will exclude 0 from the minimum search. Defaults to False.
        
        **Returns**:
        
         - The minimum value
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
        """Finds the location of the minimum value of ``param`` 
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to find minimum of. See :func:`~MARIGOLD.Condition.print_params` for options
        
        **Returns**:
        
         - ``(rstar, angle)`` location corresponding to the minimum
        """
                
        min = 10**7
        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] < min:
                    min = midas_dict[param]
                    location = (rstar, angle)

        return (location)
    
    def max_line_loc(self, param: str, angle:float) -> float:
        """Finds the location of the maximum value of ``param`` along a given line, defined by ``alpha``
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to find maximum of. See :func:`~MARIGOLD.Condition.print_params` for options
         - ``angle``: :math:`\\varphi` value to search down
        
        **Returns**:
        
         - ``rstar`` corresponding to the maximum
        """
        max = 0
        for rstar, midas_dict in self.data[angle].items():
            if midas_dict[param] > max:
                max = midas_dict[param]
                location = rstar

        return (location)
    
    def max_line(self, param: str, angle:float) -> float:
        """Finds the maximum value of ``param`` along a given line, defined by ``alpha``
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to find maximum of. See :func:`~MARIGOLD.Condition.print_params` for options
         - ``angle``: :math:`\\varphi` value to search down
        
        **Returns**:
        
         - maximum value
        """
        max = 0
        for rstar, midas_dict in self.data[angle].items():
            if midas_dict[param] > max:
                max = midas_dict[param]
                location = rstar

        return (max)
    
    def find_hstar_pos(self, method='max_dsm', void_criteria = 0.05) -> float:
        """Returns the vertical distance from the top of the pipe to the "bubble layer interface"

        .. math:: h^{*} = 1 - r^{*} \\sin{\\varphi}
        
        **Args**:
        
         - ``method``: method for determining :math:`h^{*}`, the bubble layer height. Defaults to 'max_dsm'.
             
             - ``max_dsm``, :math:`r` and :math:`\\varphi` determined by calling :any:`max_loc`\ ``('Dsm1')``. Assumes the bubble layer limit corresponds to the maximum Sauter-mean diameter
             - ``min_grad_y``, :math:`r` and :math:`\\varphi` determined by calling :any:`min_loc` for ``grad_alpha_y``, after calling :any:`calc_grad`\ ``('alpha')``
             - ``max_grad_y``:math:`r` and :math:`\\varphi` determined by calling :any:`max_loc` for ``grad_alpha_y``, after calling :any:`calc_grad`\ ``('alpha')``
             - ``max_mag_grad_y``, Not implemented
             - ``zero_void``, searches down the :math:`\\varphi = 0` line until the largest :math:`r^{*}` is found where :math:`\\alpha (r^{*})` < :math:`\\alpha_{c}`  
             - ``percent_void``, searches down the :math:`\\varphi = 0` line until the largest :math:`r^{*}` is found where :math:`\\alpha (r^{*})` < :math:`\\alpha_{c}` \\times :math:`\\alpha_{max}`  
             - ``Ryan_Ref``, :math:`h^{*} = 1 - (1.3 - 1.57 \\times 10^{-5})`, based on Ryan (2022) :math:`(r/R)_{end}`

         - ``void_criteria``: value for alpha comparisons, :math:`\\alpha_c`. Either a constant, for the ``zero_void`` option, or a percentage, for ``percent_void``. Defaults to 0.05.
        
        **Returns**:
        
         - ``h_star``, if successful. ``np.NaN`` otherwise
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
    
    def check_param(self, param:str, strict=True) -> bool:
        """method to check if a parameter is present for a condition
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter to check. See :func:`~MARIGOLD.Condition.print_params` for options
         - ``strict``: if true, will return false if a single point is missing (as long as that point isn't at r/R=1). Defaults to True.
        
        **Returns**:
        
         - ``True`` if parameter is present everywhere, ``False`` if not
        """
         
        points_to_add = []
        for angle, r_dict in self.data.items():
            param_in_angle = False
            
            for rstar, midas_dict in r_dict.items():
                if param not in midas_dict.keys():
                    if strict:
                        if abs(rstar - 1.0) < 0.001:
                            points_to_add.append( (angle, rstar))
                            pass 
                        elif points_to_add == []:
                            # return (False, rstar, angle)
                            return False
                
                else:
                    param_in_angle = True

            if (not param_in_angle) and (not strict):
                # return (False, rstar, angle)
                return False
            
            # The only points missing were at r/R = 1, which should just be 0 anyways
            if points_to_add:
                for phi, r in points_to_add:
                    self.data[phi][r][param] = 0

        return True
    
    def check_param_loc(self, param:str) -> tuple:
        """ Checks if a given parameter is present everywhere in midas_dict
        
        """

        for angle, r_dict in self.data.items():
            for rstar, midas_dict in r_dict.items():
                if param not in midas_dict.keys():
                    return (angle, rstar)

        return ()

    def calc_grad(self, param: str, recalc = False) -> None:
        """Calculates gradient of a given parameter
        
        **Args**:
        
            - ``param``: ``midas_dict`` parameter to calculate gradient of. See :func:`~MARIGOLD.Condition.print_params` for options
            - ``recalc``: recalculate when called. If False, will pull from dictionary if entry is available. Defaults to False.
        
        **Returns**:
        
            - Area averaged ``grad_'param'_total``.
            - Stores:
                - ``'grad_param_r'``
                - ``'grad_param_phi'``
                - ``'grad_param_phinor'``
                - ``'grad_param_x'``
                - ``'grad_param_y'``
                - ``'grad_param_total'``
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
                    
        return area_avg(self,grad_param_name+'_total')

    def calc_void_cov(self):
        """Calculates the void covariance
        
        .. math:: \\frac{\\langle \\alpha^2 \\rangle }{\\langle \\alpha \\rangle^2}
        
        **Returns**:
        
         - void covariance
         - stores ``self.void_cov``
        """
        # Calculates the void covariance

        # Inputs:
        #  - None

        # Stores:
        #  - "void_cov" in self

        # Returns:
        #  - void covariance

        # Mathematically performing the operation
        
        # 
        
        # 

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
        I = integrate.simpson(y=param_r_int, x=angles_int) / np.pi / area_avg(self,'alpha')**2 # Integrate wrt theta, divide by normalized area

        self.void_cov = I
        return I

    def calc_sigma_param(self, param):
        """Calculates the second moment of a parameter, :math:`\\psi`

        .. math:: \\sigma_{\\psi} = \\frac{\\langle (\\psi - \\langle \\psi \\rangle)^2 \\rangle }{\\langle \\psi \\rangle^2}
        
        **Args**: 
        
         - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        
        **Returns**:
        
         - :math:`\\sigma_{\\psi}`
        """

        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        param_avg = area_avg(self,param)

        for angle, r_dict in self.data.items():
            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                    vars_temp.append(rstar * (midas_dict[param] - param_avg)**2)
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
        I = integrate.simpson(y=param_r_int, x=angles_int) / np.pi / param_avg**2 # Integrate wrt theta, divide by normalized area
        if debug: print('Calculated sigma_alpha: ', I)

        return I

    def calc_mu3_alpha(self, param):
        """Calculates the third moment of a parameter, :math:`\\psi`
        
        .. math:: \\mu_{3, \\psi} = \\frac{\\langle (\\psi - \\langle \\psi \\rangle)^3 \\rangle }{\\langle \\psi \\rangle^3}
        
        **Args**: 
        
         - ``param``: ``midas_dict`` parameter. See :any:`print_params` for options
        
        **Returns**:
        
         - Third moment of :math:`\\psi`
        """

        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        param_avg = area_avg(self,param)

        for angle, r_dict in self.data.items():
            rs_temp = []
            vars_temp = []
            angles.append(angle * np.pi/180) # Convert degrees to radians
            for rstar, midas_dict in r_dict.items():
                if rstar >= 0:
                    rs_temp.append( rstar ) # This is proably equivalent to rs = list(r_dict.keys() ), but I'm paranoid about ordering
                    vars_temp.append(rstar * (midas_dict[param] - param_avg)**3)
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
        I = integrate.simpson(y=param_r_int, x=angles_int) / np.pi / param_avg**3 # Integrate wrt theta, divide by normalized area
        if debug: print('Calculated mu3_alpha: ', I)

        return I

    def top_bottom(self, param, even_opt='first') -> float:
        """Not sure how this isn't the same as :any:`area_avg`. 
        
        **Args**:
        
         - ``param``: ``midas_dict`` parameter. See :func:`~MARIGOLD.Condition.print_params` for options
         - ``even_opt``: for scipy.integrate.simpsons. Defaults to 'first'.
        
        **Returns**:
        
         - area-average?
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

    def calc_errors(self, param1:str, param2:str) -> float:
        """Calculates the errors, ε, between two parameters (param1 - param2) in midas_dict

        Usually want to do param1=predicted, param2=experimental

        If param2 = 0, relative errors are considered 0
        
        **Args:**

         - ``param1``: parameter to calculate error between (predicted)
         - ``param2``: parameter to calculate error between (experimental)

        Stores:
         - error, ``'eps_param1_param2'``, param1 - param2
         - relative, ``'eps_rel_param1_param2'``, (param1 - param2) / param2
         - absolute relative, ``'eps_abs_rel_param1_param2'``, | param1 - param2 | / param2
         - square, ``'eps_sq_param1_param2'``, (param1 - param2)**2
         - relative square, ``'eps_rel_sq_param1_param2'``, ((param1 - param2)/param2)**2
        
        Returns:
         - Area-averaged error
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

        return area_avg(self,param_error_name)
    
    def calc_AA_error(self, param1:str, param2:str) -> float:
        """Calculates the error, ε, between the area-average of two parameters (⟨param1⟩ - ⟨param2⟩) in ``midas_dict``
        
        **Args:**

         - ``param1``: parameter to calculate error between (predicted)
         - ``param2``: parameter to calculate error between (experimental)
        
        Returns:
         - relative error (⟨param1⟩ - ⟨param2⟩) / ⟨param1⟩
        """

        eps = area_avg(self,param1) - area_avg(self,param2)
        rel_error = eps / area_avg(self,param1)

        return rel_error

    def calc_symmetry(self, param, sym_type = 'sym90', method = 'rmse', rel_error = False):
        """Function for checking the symmetry of a condition
        
        **Args:**

         - ``param``: ``midas_dict`` parameter to plot. See :func:`~MARIGOLD.Condition.print_params` for options
         - ``sym_type``: type of symmetry to compare. Defaults to 'sym90'. Currently the only option.
         - ``method``: error options. Defaults to 'rmse'. Options:
             - ``'rmse'``, :math:`\\epsilon_{sym} = \\sqrt{\\bar{\\epsilon(r)^{2}}}`
             - ``'mean'``, :math:`\\epsilon_{sym} = \\bar{\\epsilon(r)}`
         - ``rel_error``: will divide by param. Defaults to False.
        
        Returns:
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
        """Calculate the symmetry error

        Not sure this is necessary, could use :meth:`~MARIGOLD.Condition.Condition.calc_symmetry` and :meth:`~MARIGOLD.Condition.Condition.area_avg` to achieve a similar thing
        
        99% sure this was AI generated by quanz

        **Args:**

         - ``param``: ``midas_dict`` parameter to plot. See :func:`~MARIGOLD.Condition.print_params` for options
         - ``sym_type``: _description_. Defaults to 'sym_half'.
         - ``rel_error``: _description_. Defaults to True.
         - ``even_opt``: _description_. Defaults to 'first'.
        
        Raises:
         - ``ValueError``: _description_
        
        Returns:
         - _description_
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

        return area_avg(self,'lambda')
    
def print_tab_keys() -> None:
    """Alias of print_params()
    """
    print_params()
    return

def print_params() -> None:
    """Convenience function in case you forget a parameter
     
    Generally, the available the available parameters are:

     - ``'roverR'``, nondimensional radial location
     - ``'time'``, measurement time
     - ``'frequency'``, data acquisition frequency
     - ``'num_spherical'``, number of spherical bubbles observed
     - ``'num_distorted'``, number of distorted bubbles observed
     - ``'num_cap'``, number of cap bubbles observed
     - ``'num_slug'``, number of slug bubbles observed
     - ``'num_G1'``, number of Group I bubbles observed
     - ``'num_G2'``, number of Group II bubbles observed
     - ``'num_total'``, total number of bubbles observed
     - ``'obs_0'``, number of bubbles observed at sensor 0
     - ``'obs_1'``, number of bubbles observed at sensor 1
     - ``'obs_2'``, number of bubbles observed at sensor 2
     - ``'obs_3'``, number of bubbles observed at sensor 3
     - ``'bub_freq'``, bubble freqency
     - ``'pair_spherical'``, number of spherical bubble signals paired
     - ``'pair_distorted'``, number of distorted bubble signals paired
     - ``'pair_cap'``, number of cap bubble signals paired
     - ``'pair_slug'``, number of slug bubble signals paired
     - ``'total_paired'``, total number bubble signals paired
     - ``'percent_paired'``, percentage of observed signals that were paired
     - ``'alpha_spherical'``, void fraction contribution from spherical bubbles
     - ``'alpha_distorted'``, void fraction contribution from distorted bubbles
     - ``'alpha_cap'``, void fraction contribution from cap bubbles
     - ``'alpha_slug'``, void fraction contribution from slug bubbles
     - ``'alpha_G1'``, void fraction contribution from Group I bubbles
     - ``'alpha_G2'``, void fraction contribution from Group II bubbles
     - ``'alpha'``, total void fraction 
     - ``'ai_spherical'``, interfacial area concentration contribution from spherical bubbles
     - ``'ai_distorted'``, interfacial area concentration contribution from distorted bubbles
     - ``'ai_cap'``, interfacial area concentration contribution from cap bubbles
     - ``'ai_slug'``, interfacial area concentration contribution from slug bubbles
     - ``'ai_G1'``, interfacial area concentration contribution from Group I bubbles
     - ``'ai_G2'``, interfacial area concentration contribution from Group II bubbles
     - ``'ai'``, total interfacial area concentration :math:`[m^{-1}]`
     - ``'Dsm1'``, Sauter-mean diameter from Group I bubbles :math:`[mm]`
     - ``'Lcl1'``, chord length from Group I bubbles :math:`[mm]`
     - ``'Dsm2'``, Sauter-mean diameter from Group II bubbles :math:`[mm]`
     - ``'Lcl2'``, chord length from Group II bubbles :math:`[mm]`
     - ``'ug1'``, Average interface velocity from Group I bubbles :math:`[m/s]`
     - ``'ug2'``, Average interface velocity from Group II bubbles :math:`[m/s]`
     - ``'sigma_ug1'``, standard deviation of Group I bubble velocities (?) :math:`[m/s]`
     - ``'sigma_ug2'``, standard deviation of Group II bubble velocities (?) :math:`[m/s]`
     - ``'fluctuation'`` :math:`[m/s]`
     - ``'alpha_ug1'``, Local :math:`j_{g}` calculated by :math:`\\alpha` and Group I bubble velocity :math:`[m/s]`
     - ``'alpha_ug2'``, Local :math:`j_{g}` calculated by :math:`\\alpha` and Group II bubble velocity :math:`[m/s]`
     - ``'alpha_ug'``, Total local :math:`j_{g}` :math:`[m/s]`
     - ``'alpha_Dsm1'``, Product of alpha and :math:`D_{sm,1}` :math:`[mm]`
     - ``'alpha_Dsm2'``, Product of alpha and :math:`D_{sm,2}` :math:`[mm]`
     - ``'r01'``, distance between sensor 0 and 1 :math:`[mm]`
     - ``'r02'``, distance between sensor 0 and 2 :math:`[mm]`
     - ``'r03'``, distance between sensor 0 and 3 :math:`[mm]`
     - ``'r12'``, distance between sensor 1 and 2 :math:`[mm]`
     - ``'r13'``, distance between sensor 1 and 3 :math:`[mm]`
     - ``'r23'``, distance between sensor 2 and 3 :math:`[mm]`
     - ``'vf'``, liquid velocity :math:`[m/s]`
     - ``'jf_loc'``, local superficial liquid velocity :math:`[m/s]`
     - ``'jf'``, superficial liquid veloicty :math:`[m/s]`
     - ``'delta_p'``, pressure difference measured by Pitot-static probe :math:`[psi]`
     - ``'sigma_delta_p'``, standard deviation of pressure difference measured by Pitot-static probe :math:`[psi]`
     - ``'vr'``, relative velocity :math:`[m/s]`
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
