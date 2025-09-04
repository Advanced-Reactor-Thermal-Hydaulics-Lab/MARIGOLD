from .config import *
from .Condition import *

"""
Helpful utilities for analyzing data, writing to a different format, etc. 

"""

def comp_cond(cond1:Condition, cond2:Condition, tag = 'run_ID', rmesh_preference = '1') -> Condition:
    """ Collate data from cond1 and cond2 into a single condition
    
    Each param will be tagged with "tag", options are
    - run_ID, use cond.run_ID
    - jf, use cond.jf
    - jgloc, use cond.jgloc
    - port, use cond.port
    - name, use cond.name
    - exp_cfd, tag1 -> exp, tag2 -> CFD
    - If given a tuple, tag1 -> tuple[0], tag2 -> tuple[1]

    rmesh_preference controls which mesh to use. If the mesh point doesn't exist for the other condition,
    it will be linearly interpolated.

    Options:
     - '1', cond1 mesh will be used
     - '2', cond2 mesh will be used
     - 'finer', whichever condition has a finer mesh (more angles)

    """
    compCond = Condition(cond1.jgref, cond1.jgloc, cond1.jf, cond1.theta, cond1.port, cond1.database)

    if type(tag) == tuple:
        tag1 = tag[0]
        tag2 = tag[1]

    elif tag.lower() == 'run_ID':
        tag1 = cond1.run_ID
        tag2 = cond2.run_ID

    elif tag.lower() == 'jf':
        tag1 = cond1.jf
        tag2 = cond2.jf

    elif tag.lower() == 'jgloc':
        tag1 = cond1.jf
        tag2 = cond2.jf

    elif tag.lower() == 'port':
        tag1 = cond1.port
        tag2 = cond2.port

    elif tag.lower() == 'name':
        tag1 = cond1.name
        tag2 = cond2.name

    elif tag.lower() == 'exp_cfd':
        tag1 = 'exp'
        tag2 = 'CFD'
    
    else:
        print("Invalid tag selected, defaulting to 1 and 2")
        tag1 = '1'
        tag2 = '2'

    # Determine which condition has more angles, use that one for rmesh as well
    if rmesh_preference.lower() == 'finer':
        if len(cond1.data.keys()) > len(cond1.data.keys()):
            rmesh_cond = cond1
            not_rmesh_cond = cond2
        else:
            rmesh_cond = cond2
            not_rmesh_cond = cond1

    elif rmesh_preference == '1':
        rmesh_cond = cond1
        not_rmesh_cond = cond2

    elif rmesh_preference == '2':
        rmesh_cond = cond2
        not_rmesh_cond = cond1

    compCond._angles = rmesh_cond.data.keys()
    
    for angle in compCond._angles:
        compCond.data.update({angle:{}})

        for rstar, data_dict in rmesh_cond.data[angle].items():
            compCond.data[angle].update({rstar:{}})
            for param, val1 in data_dict.items():
                label1 = param + '_' + tag1
                label2 = param + '_' + tag2

                compCond.data[angle][rstar].update({label1:val1})
                
                try:
                    val2 = not_rmesh_cond(angle*np.pi/180, rstar, param, interp_method='linear')
                    compCond.data[angle][rstar].update({label2:val2})
                except KeyError:
                    
                    compCond.data[angle][rstar].update({label2:np.nan})

    return compCond

def listdir_abs(directory):
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            yield os.path.abspath(os.path.join(dirpath, f))

def write_excel(cond):
    """A function to write an Excel sheet based on a condition object

    :param cond: Condition to write Excel sheet for
    :type cond: MARIGOLD.Condition
    """


    return

def write_pdf(cond:Condition, output_tex = None):
    """ Export data from a condition to a pdf, using LaTeX table

    calls pdflatex, so that must be installed for this to work properly
    
    """
    # try:
    #     run("pdflatex")
    # except:
    #     print("Error running pdflatex. Is it installed?")
    #     return -1
    
    if output_tex is None:
        output_tex = "temp.tex"

    with open(output_tex, 'w') as f:
        print("\
\\documentclass{article}\n \
\\nonstopmode\n\
\\usepackage{array}\n \
\\usepackage{gensymb}\n\
\\begin{document}\n \
\\begin{center}\n\
\\section*{$j_{f}$ = %0.1f, $j_{g}$ = %0.3f, Port %s} " % (cond.jf, cond.jgref, cond.port), file = f)
        for angle, r_dict in cond.data.items():
            print( "\\begin{tabular}{|m{1.5cm}|m{1.5cm}|m{1.5cm}|m{1.5cm}|m{1.5cm}|m{1.5cm}|}", file = f)
            print( " \\hline", file = f)
            print( " $\\varphi [\\degree]$ & r/R & $\\alpha [-]$ & $a_{i} [m^{-1}]$ & $v_{g} [m/s]$ & $D_{sm} [mm]$ \\\\", file = f)
            print( " \\hline\\hline", file = f)
            
            for rstar, midas_dict in r_dict.items():
                
                print(f"{angle:1.1f} & {rstar:1.1f} & {midas_dict['alpha']:0.3f} & {midas_dict['ai']:0.1f} & {midas_dict['ug1']:0.2f} & {midas_dict['Dsm1']:0.2f} \\\\", file = f)
            
            print( " \\hline\\hline", file = f)
            print("\\end{tabular}", file=f)
        
        print("\
\\end{center}\n \
\\end{document}" , file = f)

    run(f"pdflatex {output_tex} -interaction=nonstopmode")
    return 1

def write_csv(cond:Condition, output_name = None, param_list = ['alpha', 'ai', 'ug1', 'Dsm1']):
    if output_name is None:
        output_name = cond.name + ".csv"
    
    with open(output_name, 'w') as f:
        print('angle, rstar', end = '', file = f)
        for param in param_list:
            print(', ', end='', file=f)
            print(param, end = '', file = f)
        print('', file = f)
        for angle, r_dict in cond.data.items():
            for rstar, midas_dict in r_dict.items(): 
                print(f"{angle}, {rstar}", file = f, end = '')

                for param in param_list:
                    print(', ', end='', file=f)
                    print(f"{midas_dict[param]:0.3f}", end = '', file = f)
                
                print('', file=f)


    return 1
