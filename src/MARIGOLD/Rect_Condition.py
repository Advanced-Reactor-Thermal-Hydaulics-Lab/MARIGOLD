from .config import *
from.Condition import Condition
import warnings

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

        # Data is stored in this phi array. 3 layers of dictionary
        # phi [angle] gives a dictionary with the various r/R
        # phi [angle][r/R] gives a dictionary with the MIDAS output
        # So phi[angle][r/R]['alpha'] should give you the void fraction at r/R for phi = angle
        # This structure is initialized with zeros for the MIDAS output at the pipe center and wall
        #self.phi = deepcopy(dict( zip(angles, deepcopy([ {0.0: dict( zip(tab_keys, [0]*len(tab_keys)) ), 1.0: dict(zip(tab_keys, [0]*len(tab_keys)) ) } ]) * len(angles)) ))
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