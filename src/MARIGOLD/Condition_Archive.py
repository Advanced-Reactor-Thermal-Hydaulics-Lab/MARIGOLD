from .config import *
from .Condition import Condition
import warnings

class Iskandrani_Condition(Condition):
    
    def __init__(self, jf, jg):

        self.name = f"Iskandrani_jf={jf}_jg={jg}"

        self.jf = jf
        self.jg = jg

        self.data = {}

    def __eq__(self, __o: object) -> bool:
        if isinstance(__o, Iskandrani_Condition):

            return ((self.jf == __o.jf) and (self.jg == __o.jg) )
        
        return False

    def __hash__(self) -> int:
        return hash(repr(self))

    def __repr__(self) -> str:
        return self.name
    
    def area_avg(self, param):
        print("Warning, this is a shitty area-average because we don't have the whole pipe cross section")
        # Using Bottin's scheme
        integral = 0
        # Need rs at midpoints, only on the positive side
        _rs = []
        for i in range(0, len(self.data['rs'])-1):
            if self.data['rs'][i+1] < 0: break
            r_to_append = 0.5*(self.data['rs'][i+1] + self.data['rs'][i]) 
        
            _rs.append(r_to_append)
        
        #print(self.data['rs'])
        #print(_rs)    
        # Calculate areas
        _As = []
        for i, r in enumerate(_rs):
            if i == 0:
                _As.append( 1 / np.pi * (np.arccos(r) - r*np.sqrt(1 - r**2)) )
            else:
                _As.append( 1 / np.pi * (np.arccos(r) - r*np.sqrt(1 - r**2) - np.arccos(_rs[i-1]) + _rs[i-1]*np.sqrt(1 - _rs[i-1]**2) ))

        temp = _As.copy()

        temp.reverse()

        _As.append( 2* (1 / np.pi * (np.arccos(0) - 0*np.sqrt(1 - 0**2)) - 1 / np.pi * (np.arccos(_rs[-1]) - _rs[-1]*np.sqrt(1 - _rs[-1]**2) )))

        _As += temp

        #print(_As, sum(_As))

        if len(_As) != len(self.data[param]):
            print("Error: Area list is not the same length as data to 'integrate' ")
            print(f'{_As=}')
            print(f'{self.data[param]=}')

        for _A, val in zip(_As, self.data[param]):
            integral += (_A * val)

        return integral

class Rect_Condition(Condition):

    def __init__(self, jgref:float, jgloc:float, jf:float, theta:int, LoverD:float, width:float, depth:float, database:str, fluids = 'air-water', g = 9.81) -> None:
        """ Initialize Condition object

        Inputs:
        - jgref, reference superficial gas velocity
        - jgloc, local superficial gas velocity
        - jf, superficial liquid velocity
        - theta, angle of inclination of flow direction (0° is horizontal, 90° is vertical upwards)
        - LoverD, L/D location
        - width, width of channel
        - depth, depth of channel
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
        self.LoverD = LoverD
        self.width = width
        self.depth = depth
        self.database = database

        self.name = f"jf={self.jf}_jgloc={self.jgloc:0.2f}_theta={self.theta}_port={self.LoverD}_{self.database}"

        self.Dh = 4 * (width * depth) / (2*width + 2*depth)

        # Data is stored in this data array
        self.data = {}

        self.mirrored = False
        self.FR = 0 # Flow regime variable. 0 is undefined, 1 is bubbly, etc.

        self.j = self.jgloc + self.jf

        self.vwvg = -1
        self.void_cov = -1

        self.area_avg_void_sheet = -1

        # Empty dictionaries, filled when max or area avg is called
        self.area_avgs = {}
        self.circ_seg_area_avgs = {}
        self.maxs = {}
        self.mins = {}
        self._grads_calced = []

        if fluids == 'air-water':
            self.rho_f = 998       # kg/m^3
            self.rho_g = 1.204     # kg/m^3
            self.mu_f = 0.001002   # Pa s
            self.mu_g = 0.01803e-3 # Pa s

            self.sigma = 0.0728    # N/m
        
        else:
            raise NotImplementedError(f"{fluids} not available, try 'air-water'")
        
        self.g = g
        self.gz = g * np.sin(self.theta * np.pi/180)


        return

    def mirror(self, symQuad = True):
        # TODO

        pass

    def area_avg(self, param: str, even_opt='first', recalc = True) -> float:
        try:
            dummy = self.data[0][0][param]
        except KeyError as e:
            print(f"KeyError: {e}")
            if debug: print(self.data, file=debugFID)
            print(f"Could not area-average {param} for condition {self.name}")
            return None
        
        if (param in self.area_avgs.keys()) and (not recalc):
            return self.area_avgs[param] # why waste time, if we already calculated this don't do it again
        
        # We have to integrate twice, once with resepect to r, again with respect to phi
        # Start with r

        I = 0
        param_r = [] # array for parameter integrated wrt r
        xs = []
        
        if not self.mirrored:
            warnings.warn("Mirroring in area-avg")
            self.mirror()


        for x, y_dict in self.data.items():

            ys_temp = []
            vars_temp = []
            xs.append(x) # Convert degrees to radians
            for y, midas_dict in y_dict.items():
                ys_temp.append( y )
                vars_temp.append( midas_dict[param] )
            
            
            vars = [var for _, var in sorted(zip(ys_temp, vars_temp))]
            ys = sorted(ys_temp)

            if debug: print("Arrays to integrate", x, ys, vars, file=debugFID)

            if len(ys) != len(vars):
                ValueError( f"rs to integrate over {ys} must be the same length as params {vars}, occured at {x}" )
                
            try:
                param_r.append( integrate.simpson(vars, ys, even=even_opt) ) # Integrate wrt r
            except Exception as e:
                print(e)
                print(ys, vars)
            if debug: print("calculated integral:", integrate.simpson(vars, ys, even=even_opt), file=debugFID)
                #I = 2 * np.pi
        if debug: print("Integrated wrt y", param_r, file=debugFID)

        param_r = [param for _, param in sorted(zip(xs, param_r))]
        xs = sorted(xs)

        I = integrate.simpson(param_r, xs, even=even_opt) / (self.width * self.depth) # Integrate wrt theta, divide by normalized area
        self.area_avgs.update({param: I})
        return I

        pass

class Yang_Condition(Condition):
    def __init__(self, jf, jg):
        self.jf = jf
        self.jg = jg

        self.name = f'Yang_{jf}_jg'

        self.phi = {}
        self.mirrored = False
        self.theta = 0

        self._angles = np.arange(0, 361, 45) # Hardcoded 45° increments

    def __eq__(self, __o: object) -> bool:
        if isinstance(__o, Yang_Condition):

            return ((self.jf == __o.jf) and (self.jgref == __o.jg))
        
        return False

    def pretty_print(self) -> None:
        print(self.name)
        for angle, r_dict in self.phi.items():
            print(angle)
            for r, midas_output in r_dict.items():
                print(f"\t{r}")
                print("\t\t", midas_output) 
        return
