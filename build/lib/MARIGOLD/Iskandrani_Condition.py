from .config import *

class Iskandrani_Condition:
    
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
