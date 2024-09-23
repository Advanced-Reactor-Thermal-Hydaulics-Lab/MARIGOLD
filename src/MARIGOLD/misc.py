from .Condition import Condition
from .config import *
from subprocess import run
from shutil import copy2

def comp_cond(cond1:Condition, cond2:Condition, tag = 'run_ID') -> Condition:
    """ Collate data from cond1 and cond2 into a single condition
    
    Each param will be tagged with "tag", options are
    - run_ID, use cond.run_ID
    - jf, use cond.jf
    - jgloc, use cond.jgloc
    - port, use cond.port
    - name, use cond.name
    - exp_cfd, tag1 -> exp, tag2 -> CFD
    - If given a tuple, tag1 -> tuple[0], tag2 -> tuple[1]

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
    
    elif type(tag) == tuple:
        tag1 = tag[0]
        tag2 = tag[1]

    # Determine which condition has more angles, use that one for rmesh as well
    if len(cond1.data.keys()) > len(cond1.data.keys()):
        compCond._angles = cond1.data.keys()
        rmesh_cond = cond1
        not_rmesh_cond = cond2
    else:
        compCond._angles = cond2.data.keys()
        rmesh_cond = cond2
        not_rmesh_cond = cond1

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
                    
                    compCond.data[angle][rstar].update({label2:np.NaN})

    return compCond

def write_excel(cond):
    """ Export data from a condition to an excel sheet
    
    # TODO
    """


    return


def write_inp(roverR, filename, probe_number = 'AM4-5', r01=1.408, r02=1.593, r03=1.597, r12=0.570, r13=0.755, r23=0.343, directory = os.getcwd(), detailedOutput=0, signalOutput=0):
    """ Write an .inp file for MIDAS
    
    """
    with open(os.path.join(directory, 'Input.inp'), 'w') as f:
        f.write(
f"*PROBE NUMBER: {probe_number}\n\
probetype=4\n\
r/R={roverR}\n\
r01={r01}\n\
r02={r02}\n\
r03={r03}\n\
r12={r12}\n\
r13={r13}\n\
r23={r23}\n\
frequency=50000\n\
measuretime=30\n\
Filename={filename}\n\
DetailedOutput={detailedOutput}\n\
SignalOutput={signalOutput}\n\
CLhistmax=1.0\n\
CLhistbins=10000\n\
*\n\
*New Squaring Method-AMFLSquareLoc\n\
*0-original method \n\
*1-stepping backwards and checking median filtered slope\n\
*2-three criteria slope correction\n\
AMFLSquareLoc=2\n\
SqLocThresh=0.0\n\
*\n\
*BatchNum=2\n\
*\n\
*r/R=0.0\n\
*measuretime=60\n\
*Filename=r0i\n\
*\n\
*r/R=0.0\n\
*measuretime=60\n\
*Filename=r0o\n\
*\n\
end"
        )

    return

def process_dir(target_dir:str, probe_number:str, r01:float, r02:float, r03:float, r12:float, r13:float, r23:float, signal_output=0, detailed_output=0):
    """ Runs MIDAS for every dat file in a given directory

    Makes a new folder, auto_reprocessed_data_TIMESTAMP, where the .tab files will be put.

    Inputs:
     - target_dir, directory to process
     - Probe number, for identification
     - Probe measurements (r01, r02, etc.). In mm. Same as in .inp file
     - signal_output, makes the _MedianSig, _NormSig, etc. files

    Outputs:
     - Returns name of directory the reprocessed files are in
    
    """
    os.chdir(target_dir)
    current_time = datetime.now()
    timestamp = f"{current_time.month}-{current_time.day}-{current_time.year}_{current_time.hour}-{current_time.minute}"

    reprocessed_dir = os.path.join(target_dir, 'auto_reprocessed_data_'+ timestamp)
    
    os.makedirs( reprocessed_dir )

    copy2(os.path.join(target_dir, "MIDASv1.14d.exe"), reprocessed_dir)

    os.chdir(reprocessed_dir)

    for file in os.listdir(target_dir):
        # print(file)
        if file.split('.')[-1] == 'dat':
            copy2(os.path.join(target_dir, file), reprocessed_dir)

            roverR = 0.1 * int(file[1])
            write_inp(roverR, file.split('.')[0], probe_number = probe_number, r01=r01, r02=r02, r03=r03, r12=r12, r13=r13, r23=r23, directory=reprocessed_dir, signalOutput=signal_output, detailedOutput=detailed_output)

            comp_process = run(os.path.join(reprocessed_dir, 'MIDASv1.14d.exe'), cwd = reprocessed_dir, shell=True)

            if comp_process.returncode != 0:
                print(comp_process)
            try:
                os.remove(os.path.join(reprocessed_dir, file))
            except OSError as e:
                print("Failed to remove .dat file, ", e)

    return reprocessed_dir
            