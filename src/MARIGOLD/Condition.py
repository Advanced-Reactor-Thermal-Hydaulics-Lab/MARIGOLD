"""
This is the the Condition class, which is the main abstraction to handle local probe data

Data is stored in the Condition.phi property. It's actually 3 layers of dictionary
phi [angle] gives a dictionary with the various r/R
phi [angle][r/R] gives a dictionary with the MIDAS output
The MIDAS output is itself a dictionary, with the keys listed in the "tab_keys" array
So phi[angle][r/R]['alpha'] should give you the void fraction at r/R for phi = angle
This structure is initialized with zeros for the MIDAS output at the pipe center and wall

Methods:
  pretty_print- Prints out the data in the condition in a more human-readable way
  mirror- Copies any data in the negative r/R for a given phi to the corresponding complementary angle. 
          Also copies data assuming some kind of symmetry. Ensures all angles (22.5° increments) are 
          represented with data, and that the data ranges from r/R 0-1. Data is guaranteed to exist for
          at least r/R = 0 and r/R = 1 (filled with zero_data if no data exists) for plotting
  approx_vf- calculates approximate vf based on a simple power law profile
  approx_vf_Kong- calculates approximate vf based on Kong's asymmetric method (TODO)
  calc_vr- calculates the relative velocity vg - vf
  calc_vgj- calculate local vg - j
  calc_grad- calculates gradient and saves the local information in self.phi[angle][r/R]['grad_"param name"_"direction"]
             where direction can be "r", "phi", "total" or "y" as of now
  area_avg- 
  line_avg-
  line_dev-
  void_area_avg-
  calc_void_cov-
  calc_sigma_alpha-
  calc_mu3_alpha-
  top_bottom-
  plot_profiles- the 2D line plots we make 
  plot_contours- cool contour plots
  plot_surface- rad surface plots
  rough_FR_ID- rough flow regime identification
  TD_FR_ID- Flow regime identification from Taitel and Dukler (TODO)
"""

from .config import *
from scipy import interpolate

class Condition:
    debugFID = None
    def __init__(self, jgP3:float, jgloc:float, jf:float, theta:int, port:str, database:str) -> None:
        
        self.jgp3 = jgP3
        self.jf = jf
        self.jgloc = jgloc
        self.theta = theta
        self.port = port
        self.database = database

        self.name = f"jf={self.jf}_jgloc={self.jgp3}_theta={self.theta}_port={self.port}_{self.database}"

        # Data is stored in this phi array. 3 layers of dictionary
        # phi [angle] gives a dictionary with the various r/R
        # phi [angle][r/R] gives a dictionary with the MIDAS output
        # So phi[angle][r/R]['alpha'] should give you the void fraction at r/R for phi = angle
        # This structure is initialized with zeros for the MIDAS output at the pipe center and wall
        self._angles = np.arange(0, 361, 22.5) # HARDCODED 22.5 degree increments
        #self.phi = deepcopy(dict( zip(angles, deepcopy([ {0.0: dict( zip(tab_keys, [0]*len(tab_keys)) ), 1.0: dict(zip(tab_keys, [0]*len(tab_keys)) ) } ]) * len(angles)) ))
        self.phi = {}

        self.mirrored = False
        self.FR = 0 # Flow regime variable. 0 is undefined, 1 is bubbly, etc.

        self.j = self.jgloc + self.jf

        if 'D' in self.port:
            self.LoverD = int(self.port.strip('D'))
        else: # Assume it's PITA
            if self.port == 'P1':
                self.LoverD = 30
            elif self.port == 'P2':
                self.LoverD = 66
            elif self.port == 'P3':
                self.LoverD = 110
            else:
                self.LoverD = -1
                print(f"Warning: Could not determine port L/D for {self}")

        self.vwvg = -1
        self.void_cov = -1

        self.area_avg_void_sheet = -1

        if database == 'Ryan':
            self.Dh = 0.0254 # m, for Ryan
            self.marker_type = 'o'
            self.marker_color = 'r'

        elif database == 'Kong':
            self.Dh = 0.1016
            self.marker_type = 's'
            self.marker_color = 'b'

        elif database == 'Talley':
            self.Dh = 0.0381
            self.marker_type = '^'
            self.marker_color = 'g'

        elif database == 'Yadav':
            self.Dh = 0.0508
            self.marker_type = 'D'
            self.marker_color = 'purple'

        else:
            print(f"Warning: Could not determine Dh for {cond.name}")
            self.Dh = np.NaN
            self.marker_type = '$?$'
            self.marker_color = 'yellow'

        # Empty dictionaries, filled when max or area avg is called
        self.area_avgs = {}
        self.circ_seg_area_avgs = {}
        self.maxs = {}
        self.mins = {}

    def __eq__(self, __o: object) -> bool:
        if isinstance(__o, Condition):

            return ((self.jf == __o.jf) and (self.jgp3 == __o.jgp3) and (self.theta == __o.theta) and (self.port == __o.port) and self.database == __o.database)
        
        return False

    def __hash__(self) -> int:
        return hash(repr(self))

    def __repr__(self) -> str:
        return self.name

    def pretty_print(self, print_to_file= True, FID=debugFID, mirror=False) -> None:
        print(f"jf = {self.jf}\tjg = {self.jgp3}\ttheta = {self.theta}\t{self.port}\t{self.database}", file=FID)
        
        if mirror:
            self.mirror()
        
        if print_to_file:
            for angle, r_dict in self.phi.items():
                print(angle, file=FID)
                for r, midas_output in r_dict.items():
                    print(f"\t{r}", file=FID)
                    print("\t\t", midas_output, file=FID)
        else:
            for angle, r_dict in self.phi.items():
                print(angle)
                for r, midas_output in r_dict.items():
                    print(f"\t{r}")
                    print("\t\t", midas_output)
        return

    def mirror(self, sym90 = True, axisym = False) -> None:
        # Mirror data, so we have data for every angle
        #
        # Quadrant definitions:
        #
        #             phi =  90
        #              , - ~ ~ ~ - ,
        #          , '       |        ' ,
        #        ,           |            ,
        #       ,     II     |    I        ,
        #      ,             |             ,
        #  180 ,-------------|-------------, 0
        #      ,             |             ,
        #       ,    III     |   IV       ,
        #        ,           |           ,
        #          ,         |        , '
        #            ' - , _ _ _ ,  '
        #                   270

        # Only ever call this function once
        if self.mirrored:
            return

        # First step is to find the phi angles that have data        
        angles_with_data = []

        for angle, rdict in self.phi.items():
            for rstar, midas_data in rdict.items():
                if any(midas_data.values()):
                    angles_with_data.append(angle)
                    break
        if debug: print('Angles with data: ', angles_with_data, file=debugFID)

        # Next, take any negative values for the existing data and copy them to the 
        # +180 phi angle (the complementary angle)
        # Also delete the negative data from the existing phi, so all angles only have postive entries
        angles_to_add = []
        for angle in angles_with_data:
            if angle <= 180:
                comp_angle = angle + 180

            if (comp_angle not in angles_with_data) and (comp_angle <= 360) and (comp_angle not in angles_to_add):

                data = deepcopy(self.phi[angle])
                rs = list(deepcopy(self.phi[angle]).keys())

                for r in rs:
                    if r > 0:
                        data.pop(r)
                
                for r in rs:
                    if r < 0:
                        data[-r] = self.phi[angle].pop(r)
                        data.pop(r)
                    
                    elif r == 0:
                        pass
                
                self.phi[angle].update({1.0: zero_data}) 

                # There should always be data at r/R 0 so we can plot contours
                try: 
                    dummy = self.phi[angle][0.0]
                except:
                    self.phi[angle].update({0.0: zero_data})  # just in case


                self.phi.update({comp_angle: {}})
                self.phi[comp_angle].update( {1.0: zero_data} )
                self.phi[comp_angle].update( data )
                

                angles_to_add.append(comp_angle)

        angles_with_data += angles_to_add
        if debug: print('Angles with data after comp_angle: ', angles_with_data, file=debugFID)

        if (360 not in (angles_with_data)) and (0 in angles_with_data):
            ref_angle = 0
            data = deepcopy(self.phi[ref_angle])
            self.phi.update({360: {}})
            self.phi[360].update( data )
            angles_with_data.append(360)

        # Now comes the actual mirroring step. Need data for every angle, incremements of 22.5° (self._angles)
        if (self.theta == 90 or axisym) and (self.port != 'P4' and self.port != 'P5'): 
            # axisymmetric
            for angle in self._angles:
                if angle not in angles_with_data:
                    ref_angle = angles_with_data[0]
                    data = deepcopy( self.phi[ref_angle] )
                
                data.update({1.0: zero_data}) # Cuz it'll get popped in the second loop. Also paranoia
                self.phi.update({angle: {}})
                self.phi[angle].update( data )

        elif sym90: 
            # symmetric across the 90 degree line
            for angle in self._angles:
                if angle not in angles_with_data:

                    data = {0.0: zero_data, 1.0: zero_data} # Fine if this gets overwritten, just need to make sure there's some data at 0 for plotting

                    if angle <= 90:
                        # Quadrant I, should usually have data here, but if we don't, try to copy data from Q2
                        ref_angle = 180 - angle
                        try:
                            data = deepcopy(self.phi[ref_angle])
                        except KeyError:
                            if debug: print(f"No data found for {angle} when mirroring {self.name}, defaulting to 0s")
                            data = {0.0: zero_data}

                    elif angle > 90 and angle <= 180:
                       # Quadrant II, mirror from Quadrant I
                        ref_angle = 180 - angle
                        data = deepcopy(self.phi[ref_angle])

                    elif angle > 180 and angle <= 270:
                        # Quadrant III, should be covered by the negative of Quadrant I
                        # But if we're here there's no data here. So make sure it's 0
                        ref_angle = 540 - angle
                        try:
                            data = deepcopy(self.phi[ref_angle])
                        except KeyError:
                            if debug: print(f"No data found for {angle} when mirroring {self.name}, defaulting to 0s")
                            data = {0.0: zero_data}

                    elif angle > 270 and angle < 360:
                        # Quadrant IV, mirror from Quadrant III
                        ref_angle = 540 - angle
                        data = deepcopy(self.phi[ref_angle])

                    elif angle == 360:
                        ref_angle = 0
                        data = deepcopy(self.phi[ref_angle])

                    data.update({1.0: zero_data}) # paranoia
                    
                    # Check if data exists at zero, and if not, just put some zero data in
                    try: 
                        dummy = data[0.0]
                    except:
                        data.update({0.0: zero_data}) # just in case
                    
                    if angle > 360: continue # Just in case
                    self.phi.update({angle: {}})
                    self.phi[angle].update( data )

            # Check if data exists at zero, and if not, just put some zero data in
                try: 
                    dummy = self.phi[angle][0.0]
                except:
                    self.phi[angle].update({0.0: zero_data}) # just in case


        else:
            # No symmetry being assumed. But we still want data at every angle. If it doesn't exist, must be 0
            for angle in self._angles:
                if angle not in angles_with_data:
                    data = {0.0: zero_data}
                    data.update({1.0: zero_data})
                    if angle > 360: continue
                    self.phi.update({angle: {}})
                    self.phi[angle].update( data )

        self.mirrored = True
        return
    
    def approx_vf(self, n=7):
        # If necessary, approximate vf using power law relation
        self.mirror()

        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    dummy = midas_dict['vf']
                    if debug: print(f"approx_vf: data found for {angle}\t{rstar}", file=debugFID)
                except:
                    vf_approx = (n+1)*(2*n+1) / (2*n*n) * (self.jf / (1-self.area_avg('alpha'))) * (1 - abs(rstar))**(1/n)
                    midas_dict.update({'vf': vf_approx})

        return
    
    def approx_vf_Kong(self, n=7):
        # TODO
        self.mirror()

        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                vf_approx = (n+1)*(2*n+1) / (2*n*n) * (self.jf / (1-self.area_avg('alpha'))) * (1 - abs(rstar))**(1/n)
                midas_dict.update({'vf': vf_approx})

        return
    
    def calc_vr(self):
        self.mirror()

        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    dummy = midas_dict['vf']
                except:
                    print("Warning: Approximating vf in calculating vr, since no data found")
                    self.approx_vf()
                vg = midas_dict['ug1']
                if vg == 0: # should be the same as α = 0, could maybe switch this to that
                    vr = 0 # this is an assumption, similar to void weighting
                else:
                    vr = midas_dict['ug1'] - midas_dict['vf']
                midas_dict.update({'vr': vr})

        return

    def calc_vgj(self):
        self.mirror()

        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                try:
                    dummy = midas_dict['vf']
                except:
                    print("Warning: Approximating vf in calculating local j, since no data found")
                    self.approx_vf()
                
                j_local = midas_dict['alpha'] * midas_dict['ug1'] + (1 - midas_dict['alpha']) * midas_dict['vf']
                vgj = midas_dict['ug1'] - j_local
                midas_dict.update({'vgj': vgj})

        return

    def calc_grad(self, param: str) -> None:

        if not self.mirrored: self.mirror()

        # calculate gradient, return 
        # stores grad_params in midas dict, along with their sum (d(param)/dr, d(param)/dφ)

        grad_param_name = 'grad_' + param

        phis = self._angles.tolist()
        phis.sort()
        maxj = len(phis)

        for phi_angle, r_dict in self.phi.items():
            rs = list(r_dict.keys())

            rs.sort()
            #print(rs)

            for i in range(1, len(rs) ):
                grad_r_param = (r_dict[rs[i]][param] - r_dict[rs[i-1]][param]) / (rs[i] - rs[i-1])
                
                
                j = phis.index(phi_angle)
                try:
                    hi = self.phi[phis[(j+1) % maxj]][rs[i]][param]
                except KeyError as e:
                    if debug: print(f"Key error found when indexing {e} for hi. Likely a case of the data being zero for the adjacent point, setting to 0...", file=debugFID)
                    hi = 0

                try:
                    lo = self.phi[phis[(j-1) % maxj]][rs[i]][param]
                except KeyError as e:
                    if debug: print(f"Key error found when indexing {e} for lo. Likely a case of the data being zero for the adjacent point, setting to 0...", file=debugFID)
                    lo = 0

                if rs[i]> 0:
                    grad_phi_param = 1./rs[i] * (hi - lo) / (2* (phis[(j+1) % maxj] - phis[(j-1) % maxj]) * np.pi/180)
                else:
                    grad_phi_param = 0 # I guess?

                r_dict[rs[i]].update( {grad_param_name+'_r': grad_r_param } )
                r_dict[rs[i]].update( {grad_param_name+'_phi': grad_phi_param } )
                r_dict[rs[i]].update( {grad_param_name+'_y': grad_r_param * np.sin(phi_angle) + np.cos(phi_angle)/rs[i]*grad_phi_param } )
                r_dict[rs[i]].update( {grad_param_name+'_x': grad_r_param * np.cos(phi_angle) - np.sin(phi_angle)/rs[i]*grad_phi_param } )
                r_dict[rs[i]].update( {grad_param_name+'_total': grad_r_param+grad_phi_param } )

            # Acount for not having data at 0, average value at r/R = 0.1 and r/R = -0.1
            r_dict[0.0].update( {grad_param_name+'_r': 0.5 * grad_r_param } )
            r_dict[0.0].update( {grad_param_name+'_phi': 0.5 * grad_phi_param } )
            r_dict[0.0].update( {grad_param_name+'_y': 0.5 * (grad_r_param * np.sin(phi_angle) + np.cos(phi_angle)/rs[i]*grad_phi_param) } )
            r_dict[0.0].update( {grad_param_name+'_x': 0.5 * (grad_r_param * np.cos(phi_angle) - np.sin(phi_angle)/rs[i]*grad_phi_param) } )
            r_dict[0.0].update( {grad_param_name+'_total': 0.5 * grad_r_param+grad_phi_param } )
            
            if phi_angle <= 180:
                comp_angle = phi_angle + 180
            else:
                comp_angle = phi_angle - 180

            self.phi[comp_angle][0.0].update({grad_param_name+'_r': 0.5 * grad_r_param })
            self.phi[comp_angle][0.0].update({grad_param_name+'_phi': 0.5 * grad_phi_param })
            self.phi[comp_angle][0.0].update({grad_param_name+'_total': 0.5 * grad_r_param+grad_phi_param })

        return
    
    def fit_spline(self, param: str) -> None:
        try: dummy = self.spline_interp
        except:
            self.spline_interp = {}
        self.mirror()
        rs = []
        phis = []
        vals = []
        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items(): 
                rs.append(rstar)
                phis.append(angle * np.pi / 180)
                vals.append(midas_dict[param])
        spline_interpolant = interpolate.SmoothBivariateSpline(phis, rs, vals)
        self.spline_interp.update({param: spline_interpolant})
        return
    
    def max(self, param: str, recalc=False) -> float:
        if (param in self.maxs.keys()) and (not recalc):
            return self.maxs(param) # why waste time 
        max = 0
        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] > max:
                    max = midas_dict[param]
                    location = rstar
        self.maxs.update({param:max})
        return (max)

    def max_loc(self, param: str)-> tuple:
        max = 0
        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] > max:
                    max = midas_dict[param]
                    location = (rstar, angle)

        return (location)
    
    def min(self, param: str, recalc = False)-> float:
        if (param in self.mins.keys()) and (not recalc):
            return self.mins[param] # why waste time 
        min = 10**7
        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] < min:
                    min = midas_dict[param]
                    location = rstar
        self.mins.update({param:min})
        return (min)

    def min_loc(self, param: str)-> float:
        min = 10**7
        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                if midas_dict[param] < min:
                    min = midas_dict[param]
                    location = (rstar, angle)

        return (location)

    def min_nonzero(self, param: str)-> float:
        min = 10**7
        for angle, r_dict in self.phi.items():
            for rstar, midas_dict in r_dict.items():
                if (midas_dict[param] < min) and (midas_dict[param] != 0):
                    min = midas_dict[param]
                    location = rstar

        return (min)
    
    def max_line_loc(self, param: str, angle) -> float:
        max = 0
        for rstar, midas_dict in self.phi[angle].items():
            if midas_dict[param] > max:
                max = midas_dict[param]
                location = rstar

        return (location)
    
    def max_line(self, param: str, angle) -> float:
        max = 0
        for rstar, midas_dict in self.phi[angle].items():
            if midas_dict[param] > max:
                max = midas_dict[param]
                location = rstar

        return (max)
    
    def find_hstar_pos(self, method='max_dsm', void_criteria = 0.05) -> float:
        """ Returns the vertical distance from the top of the pipe to the bubble layer interface, as determined by the selected method.
            Void criteria = minimum void for "zero_void" mode, or % of maximum void on line for "percent_void" mode """

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
            for rstar, midas_dict in self.phi[90].items():
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
            for rstar, midas_dict in self.phi[90].items():
                if midas_dict['alpha'] < void_criteria:
                    if rstar > roverRend and rstar < max_loc:
                        roverRend = rstar
            
            if roverRend == -1: # Didn't find it on the positive side
                for rstar, midas_dict in self.phi[270].items():
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
            for rstar, midas_dict in self.phi[90].items():
                if midas_dict['alpha'] < max_void * void_criteria:
                    if rstar > roverRend and rstar < max_loc:
                        roverRend = rstar
            
            if roverRend == -1:
                for rstar, midas_dict in self.phi[180].items():
                    if midas_dict['alpha'] < max_void * void_criteria:
                        if -rstar > roverRend:
                            roverRend = -rstar

            self.roverRend = roverRend
            h_star = 1 - self.roverRend
            return h_star
        
        elif method == 'Ryan_Ref':
            try:
                dummy = self.Ref
            except:
                print(f'Condition {self.name} has no Ref value, calculating assuming water')
                self.Ref = 998 * cond.jf * cond.Dh / 0.001
            
            self.roverRend = 1.3 - 1.57e-5 * self.Ref
            h_star = 1 - self.roverRend
            return h_star

        print('Invalid method for find_h_pos')
        return np.NaN

    def area_avg(self, param: str, even_opt='first', recalc = False) -> float:
        
        # Check that the parameter that the user requested exists
        try:
            dummy = self.phi[90][1.0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}")
            return
        
        if (param in self.area_avgs.keys()) and (not recalc):
            return self.area_avgs[param] # why waste time, if we already calculated this don't do it again
        
        # We have to integrate twice, once with resepect to r, again with respect to phi
        # Start with r

        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        self.mirror()


        for angle, r_dict in self.phi.items():

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
                
            param_r.append( integrate.simpson(vars, rs, even=even_opt) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(vars, rs, even=even_opt), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)

        param_r = [param for _, param in sorted(zip(angles, param_r))]
        angles = sorted(angles)

        I = integrate.simpson(param_r, angles, even=even_opt) / np.pi # Integrate wrt theta, divide by normalized area
        self.area_avgs.update({param: I})
        return I

    def circ_segment_area_avg(self, param:str, hstar:float, ngridr=25, ngridphi=25, int_err = 10**-4) -> float:
        from scipy import interpolate
        # area averages over the circular segment with height h
        # For the bubble layer region
        # to smooth it out, integrate over an interpolated mesh
        # Check that the parameter that the user requested exists
        try:
            dummy = self.phi[90][1.0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}")
            return
        
        self.mirror()

        rs = []
        phis = []
        vals = []
        for angle, r_dict in self.phi.items():
            for r_star, midas_dict in r_dict.items():
                rs.append(r_star)
                phis.append(angle * np.pi/180) # Convert degrees to radians
                try:
                    vals.append(midas_dict[param] * r_star) # Don't forget r dr dθ
                except:
                    vals.append(np.NaN * r_star) # Don't forget r dr dθ
                    print(f"Could not find {param} for φ = {phi_angle}, r = {r_star}. Substituting NaN")

        # Set up interpolation
        rs = np.asarray(rs)
        phis = np.asarray(phis)
        vals = np.asarray(vals)   

        interp = interpolate.LinearNDInterpolator(list(zip(rs, phis)), vals)

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
        from scipy import interpolate
        # area averages over the circular segment with height h
        # For the bubble layer region
        # to smooth it out, integrate over an interpolated mesh
        # Check that the parameter that the user requested exists
        try:
            dummy = self.phi[90][1.0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}")
            return
        
        self.mirror()

        rs = []
        phis = []
        vals = []
        denom = []
        for angle, r_dict in self.phi.items():
            for r_star, midas_dict in r_dict.items():
                rs.append(r_star)
                phis.append(angle * np.pi/180) # Convert degrees to radians
                denom.append(r_star * midas_dict['alpha'])
                try:
                    vals.append(midas_dict[param] * r_star * midas_dict['alpha']) # Don't forget r dr dθ
                except:
                    vals.append(np.NaN * r_star) # Don't forget r dr dθ
                    print(f"Could not find {param} for φ = {phi_angle}, r = {r_star}. Substituting NaN")

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

        # Check that the parameter that the user requested exists
        self.mirror()

        if phi_angle not in self.phi.keys():
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}\nData for {phi_angle} not found after mirroring!")
            return


        try:
            dummy = self.phi[90][1.0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}")
            return

        r_for_int = []
        var_for_int = []

        for rstar, midas_dict in self.phi[phi_angle].items():
            if rstar not in r_for_int:
                r_for_int.append(rstar)
                var_for_int.append(midas_dict[param])

        if phi_angle <=180:
            comp_angle = phi_angle+180
            
        else:
            comp_angle = phi_angle - 180
        
        for rstar, midas_dict in self.phi[comp_angle].items():
            if rstar not in r_for_int:
                r_for_int.append(-rstar)
                var_for_int.append(midas_dict[param])

        var_for_int = [param for _, param in sorted(zip(r_for_int, var_for_int))]
        r_for_int = sorted(r_for_int)

        I = integrate.simpson(var_for_int, r_for_int, even=even_opt) / 2 # Integrate wrt theta, divide by normalized length

        return I

    def line_avg_dev(self, param:str, phi_angle:float, even_opt='first') -> float:

        # Check that the parameter that the user requested exists
        self.mirror()

        if phi_angle not in self.phi.keys():
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}\nData for {phi_angle} not found after mirroring!")
            return


        try:
            dummy = self.phi[90][1.0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}")
            return

        r_for_int = []
        var_for_int = []

        for rstar, midas_dict in self.phi[phi_angle].items():
            if rstar not in r_for_int:
                r_for_int.append(rstar)
                var_for_int.append((midas_dict[param] - self.area_avg(param))**2)

        if phi_angle <=180:
            comp_angle = phi_angle+180
            
        else:
            comp_angle = phi_angle - 180
        
        for rstar, midas_dict in self.phi[comp_angle].items():
            if rstar not in r_for_int:
                r_for_int.append(-rstar)
                var_for_int.append((midas_dict[param] - self.area_avg(param))**2)

        var_for_int = [param for _, param in sorted(zip(r_for_int, var_for_int))]
        r_for_int = sorted(r_for_int)

        I = integrate.simpson(var_for_int, r_for_int, even=even_opt) / 2 / self.area_avg(param)**2 # Integrate wrt theta, divide by normalized length

        return I


    def void_area_avg(self, param: str, even_opt='first') -> float:
        
        # Check that the parameter that the user requested exists
        try:
            dummy = self.phi[90][1.0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}")
            return
        
        # We have to integrate twice, once with resepect to r, again with respect to phi
        # Start with r

        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        self.mirror()


        for angle, r_dict in self.phi.items():

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
                
            param_r.append( integrate.simpson(vars, rs, even=even_opt) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(vars, rs, even=even_opt), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)

        param_r = [param for _, param in sorted(zip(angles, param_r))]
        angles = sorted(angles)

        I = integrate.simpson(param_r, angles, even=even_opt) / np.pi / self.area_avg('alpha') # Integrate wrt theta, divide by normalized area
        return I

    def calc_void_cov(self):

        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        for angle, r_dict in self.phi.items():
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
                
            param_r.append( integrate.simpson(vars, rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(vars, rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)
        param_r_int = [var for _, var in sorted(zip(angles, param_r))]
        angles_int = sorted(angles)
        I = integrate.simpson(param_r_int, angles_int) / np.pi / self.area_avg('alpha')**2 # Integrate wrt theta, divide by normalized area

        self.void_cov = I
        return I

    def calc_sigma_alpha(self):
        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        alpha_avg = self.area_avg('alpha')

        for angle, r_dict in self.phi.items():
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
                
            param_r.append( integrate.simpson(vars, rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(vars, rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)
        param_r_int = [var for _, var in sorted(zip(angles, param_r))]
        angles_int = sorted(angles)
        I = integrate.simpson(param_r_int, angles_int) / np.pi / alpha_avg**2 # Integrate wrt theta, divide by normalized area
        if debug: print('Calculated sigma_alpha: ', I)

        self.sigma_alpha = I
        return I
    
    def calc_mu3_alpha(self):
        I = 0
        param_r = [] # integrated wrt r
        angles = []
        
        self.mirror()

        alpha_avg = self.area_avg('alpha')

        for angle, r_dict in self.phi.items():
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
                
            param_r.append( integrate.simpson(vars, rs) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(vars, rs), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)
        param_r_int = [var for _, var in sorted(zip(angles, param_r))]
        angles_int = sorted(angles)
        I = integrate.simpson(param_r_int, angles_int) / np.pi / alpha_avg**3 # Integrate wrt theta, divide by normalized area
        if debug: print('Calculated mu3_alpha: ', I)

        self.mu3_alpha = I
        return I

    def top_bottom(self, param, even_opt='first') -> float:
        
        # Check that the parameter that the user requested exists
        try:
            dummy = self.phi[90][1.0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.phi, file=debugFID)
            print(f"Cound not area-average {param} for condition {self.name}")
            return
        
        # We have to integrate twice, once with resepect to r, again with respect to phi
        # Start with r

        I = 0
        param_r = [] # array for parameter integrated wrt r
        angles = []
        
        self.mirror()

        for angle, r_dict in self.phi.items():

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
                
            param_r.append( integrate.simpson(vars, rs, even=even_opt) ) # Integrate wrt r
            if debug: print("calculated integral:", integrate.simpson(vars, rs, even=even_opt), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt r", param_r, file=debugFID)

        param_r = [param for _, param in sorted(zip(angles, param_r))]
        angles = sorted(angles)

        I = integrate.simpson(param_r, angles, even=even_opt) / np.pi # Integrate wrt theta, divide by normalized area
        return I

    

    def calc_vwvg(self):
        print("This guy needs work, probably don't want to use it")
        self.vwvg = self.jgloc / self.area_avg('alpha')
        return
    
    def plot_profiles(self, param, save_dir = '.', show=True, x_axis='r', 
                      const_to_plot = [90, 67.5, 45, 22.5, 0], include_complement = True, 
                      rotate=False, fig_size=4, title=True) -> None:
        """Plot profiles of param over x_axis, for const_to_plot, i.e. α over r/R for φ = [90, 67.5 ... 0]. include_complement will continue with the negative side if x_axis = 'r' """
        self.mirror()
        plt.rcParams.update({'font.size': 12})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "dejavuserif"

        if rotate:

            # Set up the figure to be rotated by theta
            import matplotlib as mpl
            from matplotlib.transforms import Affine2D
            import mpl_toolkits.axisartist.floating_axes as floating_axes
            fig = plt.figure(figsize=(fig_size, fig_size))
            if x_axis == 'r':
                plot_extents = self.min(param), self.max(param)*1.1, -1, 1
                transform = Affine2D().scale(fig_size / (self.max(param)*1.1 - self.min(param)), fig_size / (1 - -1)).rotate_deg(self.theta)
            else:
                plot_extents = self.min(param), self.max(param)*1.1, 0, 360
                transform = Affine2D().scale(fig_size / (self.max(param)*1.1 - self.min(param)), fig_size / (360-0)).rotate_deg(self.theta)
            
            
            helper = floating_axes.GridHelperCurveLinear(transform, plot_extents)
            fake_ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=helper)
            ax = fake_ax.get_aux_axes(transform)

        else:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300, layout='compressed')
            fake_ax = ax

        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        # Tick marks facing in
        ax.tick_params(direction='in',which='both')

        ms = marker_cycle()
        cs = color_cycle()

        if x_axis == 'r':
            for angle in const_to_plot:
                r_dict = self.phi[angle]
                rs = []
                vals = []
                for r, midas_output in r_dict.items():
                    rs.append(r)
                    vals.append( midas_output[param])

                if include_complement:
                    if angle > 180:
                        print("Error: Cannot find complement to angle > 180. Skipping")
                    else:
                        r_dict = self.phi[angle+180]
                        for r, midas_output in r_dict.items():
                            rs.append(-r)
                            vals.append( midas_output[param])

                vals = [var for _, var in sorted(zip(rs, vals))]
                rs = sorted(rs)
                    
                ax.plot(vals, rs, label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
            ax.set_ylim(-1, 1)
            ax.set_xlim(self.min(param), self.max(param)*1.2)
            
        
        elif x_axis == 'phi':
            ax.plot([], [], label=r'$r/R$', color='white', linestyle = None)
            for rtarget in const_to_plot:
                phis = []
                vals = []
                for angle, r_dict in self.phi.items():
                    for rstar, midas_output in r_dict.items():
                        if abs(rstar - rtarget) < 0.001:
                            phis.append(angle)
                            vals.append( midas_output[param])

                vals = [var for _, var in sorted(zip(phis, vals))]
                phis = sorted(phis)
                if rotate:
                    ax.plot(vals, phis, label=f'{rtarget:0.1f}', color=next(cs), marker=next(ms), linestyle = '--')
                else:
                    ax.plot(phis, vals, label=f'{rtarget:0.1f}', color=next(cs), marker=next(ms), linestyle = '--')
            if rotate:
                ax.set_ylim(0, 360)
                ax.set_xlim(self.min(param), self.max(param))
                
            else:
                ax.set_xlim(0, 360)
                ax.set_ylim(self.min(param), self.max(param))
                
        else:
            print(f"invalid axis for plot_profiles: {x_axis}. Current supported options are 'r' and 'phi'")
            return
        
        if x_axis == 'r':
            if param == 'alpha':
                fake_ax.set_xlabel(r'$\alpha$ [-]')
            elif param == 'ai':
                fake_ax.set_xlabel(r'$a_{i}$ [1/m]')
            elif param == 'ug1':
                fake_ax.set_xlabel(r'$v_{g}$ [m/s]')
            else:
                fake_ax.set_xlabel(param)
            fake_ax.set_ylabel(r'$r/R$ [-]')
            fake_ax.set_yticks(np.arange(-1, 1.01, 0.2))
            #fake_ax.set_xticks(np.linspace(self.min(param), self.max(param), 7))

        elif x_axis == 'phi':
            if not rotate:
                if param == 'alpha':
                    fake_ax.set_ylabel(r'$\alpha$ [-]')
                elif param == 'ai':
                    fake_ax.set_ylabel(r'$a_{i}$ [1/m]')
                elif param == 'ug1':
                    fake_ax.set_ylabel(r'$v_{g}$ [m/s]')
                else:
                    fake_ax.set_ylabel(param)
                fake_ax.set_xlabel(r'$\varphi$ [-]')

                fake_ax.set_xticks([0, 90, 180, 270, 360])
                #fake_ax.set_yticks(np.linspace(self.min(param), self.max(param), 7))
            else:
                if param == 'alpha':
                    fake_ax.set_xlabel(r'$\alpha$ [-]')
                elif param == 'ai':
                    fake_ax.set_xlabel(r'$a_{i}$ [1/m]')
                elif param == 'ug1':
                    fake_ax.set_xlabel(r'$v_{g}$ [m/s]')
                else:
                    fake_ax.set_xlabel(param)
                fake_ax.set_ylabel(r'$\varphi$ [-]')
                
                fake_ax.set_yticks([0, 90, 180, 270, 360])
                # fake_ax.set_xticks(np.linspace(self.min(param), self.max(param), 7))

                #fake_ax.tick_params(axis='both', labelrotation=-self.theta)
        
        if title:
            ax.set_title(self.name)

        ax.spines['bottom'].set_position(('data', 0))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        if rotate:
            fig.add_subplot(fake_ax)
        ax.legend(loc='lower right', edgecolor='white')

        fake_ax.set_aspect('auto', adjustable='datalim', share=True)
        ax.set_aspect('auto', adjustable='datalim', share=True)
        #fake_ax.set_box_aspect(1)
        
        if show:
            plt.show()
        else:
            plt.savefig(os.path.join(save_dir, f'{param}_profile_vs_{x_axis}_{self.name}.png'))
            plt.close()
        return    

    def plot_contour(self, param:str, save_dir = '.', show=True, set_max = None, set_min = None, fig_size = 4,
                     rot_angle = 0, ngridr = 50, ngridphi = 50, colormap = 'viridis', num_levels = 100,
                     annotate_h = False, cartesian = False, h_star_kwargs = {'method': 'max_dsm', 'min_void': '0.05'}) -> None:
        
        if cartesian:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300)
        else:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300, subplot_kw=dict(projection='polar'))
        plt.rcParams.update({'font.size': 12})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"
        
        self.mirror()
        rs = []
        phis = []
        vals = []

        for phi_angle in self._angles:
            r_dict = self.phi[phi_angle]
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
        
        if cartesian:
            plt.contourf(XI, YI, parami, levels = num_levels, vmin = set_min, vmax = set_max, cmap = colormap)

            x = np.linspace(-1, 1, 100)
            plt.plot(x, np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
            plt.plot(x, -np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
        else:
            plt.contourf(PHII, RI, parami, levels = num_levels, vmin = set_min, vmax = set_max, cmap = colormap)

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

        plt.colorbar(label=param)
        plt.title(self.name)

        plt.tight_layout()
        
        #cb = plt.colorbar(ticks = [0, 0.05, 0.1, 0.15, 0.2])

        if show:
            plt.show()
        else:
            plt.savefig( os.path.join(save_dir, f"{param}_contours_{self.name}.png") )
            plt.close()
        return

    def plot_surface(self, param:str, save_dir = '.', show=True, rotate_gif=False, elev_angle = 145, 
                     azim_angle = 0, roll_angle = 180, title=True, ngridr = 50, ngridphi = 50, 
                     plot_surface_kwargs = {}) -> None:

        plt.rcParams.update({'font.size': 12})
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"

        self.mirror()
        rs = []
        phis = []
        vals = []

        for phi_angle, r_dict in self.phi.items():
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

        if 'cmap' not in plot_surface_kwargs.keys(): 
            plot_surface_kwargs.update({'cmap': 'viridis'})

        surf = ax.plot_surface(Xi, Yi, parami, **plot_surface_kwargs)
        

        #plt.legend()
        ax.set_xlabel (r'$x/R$ [-]')
        ax.set_ylabel(r'$y/R$ [-]')

        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        
        ax.set_zlim([plot_surface_kwargs['vmin'], plot_surface_kwargs['vmax']])

        if param == 'alpha':
            ax.set_zlabel(r'$\alpha$ [-]')
            fig.colorbar(surf, label= r'$\alpha$ [-]')
        else:
            ax.set_zlabel(param)
            fig.colorbar(surf, label=param)
        if title: plt.title(self.name)

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

    def rough_FR_ID(self) -> None:

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
        return

    def TD_FR_ID(self) -> None:
        #dpdxL

        return


def color_cycle():
    var_list = ['#0000FF',
                '#FF0000',
                '#00FF00',
                '#00FFFF',
                '#7F00FF',
                '#7FFF7F',
                '#007F7F',
                '#7F007F',
                '#7F7F7F',
                '#000000']
    i = 0
    while True:
        yield var_list[ i % len(var_list)]
        i += 1

def marker_cycle():
    var_list = ['o', '^', 's', 'v', 'D']
    i = 0
    while True:
        yield var_list[ i % len(var_list)]
        i += 1


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
    'jf_loc'
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