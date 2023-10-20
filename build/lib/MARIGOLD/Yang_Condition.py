from .config import *
from .Condition import Condition

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

            return ((self.jf == __o.jf) and (self.jgp3 == __o.jg))
        
        return False

    def pretty_print(self) -> None:
        print(self.name)
        for angle, r_dict in self.phi.items():
            print(angle)
            for r, midas_output in r_dict.items():
                print(f"\t{r}")
                print("\t\t", midas_output) 
        return
