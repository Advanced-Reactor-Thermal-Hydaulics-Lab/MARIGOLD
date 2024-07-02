from .Condition import Condition
from .config import *

def comp_cond(cond1:Condition, cond2:Condition, tag = 'run_ID') -> Condition:
    """ Collate data from cond1 and cond2 into a single condition
    
    Each param will be tagged with "tag", options are
    - run_ID, use cond.run_ID
    - jf, use cond.jf
    - jgloc, use cond.jgloc
    - port, use cond.port
    - name, use cond.name
    - exp_cfd, tag1 -> exp, tag2 -> CFD

    """
    compCond = Condition(cond1.jgref, cond1.jgloc, cond1.jf, cond1.theta, cond1.port, cond1.database)

    if tag == 'run_ID':
        tag1 = cond1.run_ID
        tag2 = cond2.run_ID

    elif tag == 'jf':
        tag1 = cond1.jf
        tag2 = cond2.jf

    elif tag == 'jgloc':
        tag1 = cond1.jf
        tag2 = cond2.jf

    elif tag == 'port':
        tag1 = cond1.port
        tag2 = cond2.port

    elif tag == 'name':
        tag1 = cond1.name
        tag2 = cond2.name

    elif tag == 'exp_CFD':
        tag1 = 'exp'
        tag2 = 'CFD'

    # Determine which condition has more angles, use that one for rmesh as well
    if len(cond1.data.keys()) > len(cond1.data.keys()):
        compCond._angles = cond1.data.keys()
        rmesh_cond = cond1
    else:
        compCond._angles = cond2.data.keys()
        rmesh_cond = cond2

    for angle in compCond._angles:
        compCond.data.update({angle:{}})

        for rstar, data_dict in rmesh_cond[angle].items():
            compCond[angle].update({rstar:{}})
            for param, val1 in data_dict.items():
                label1 = param + '_' + tag1
                label2 = param + '_' + tag2

                compCond[angle][rstar].update({label1:val1})
                
                val2 = cond2(angle*np.pi/180, rstar, param, interp_method='linear')
                compCond[angle][rstar].update({label2:val2})

    return compCond

def write_excel(cond):
    """ Export data from a condition to an excel sheet
    
    # TODO
    """


    return