from .Condition import Condition
from .config import *
from .operations import *
from subprocess import run
from shutil import copy2

"""
This file contains functions related to interacting with raw data files from the NI card, MIDAS or PP

- write_inp, for MIDAS
- write_pitot_inp, for PP
- tdms_to_dat, converts an NI .tdms file to a .dat file
- process_dir, processes entire directory

"""


def write_inp(roverR, filename, probe_number = 'AM4-5', r01=1.408, r02=1.593, r03=1.597, r12=0.570, r13=0.755, r23=0.343, directory = os.getcwd(), detailedOutput=0, signalOutput=0, inp_name = 'Input.inp', measure_time = 30):
    """ Write an .inp file for MIDAS
    
    """
    with open(os.path.join(directory, inp_name), 'w') as f:
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
measuretime={measure_time}\n\
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

def write_pitot_inp(roverR, filename, URV = 5.248, LRV = 1.054, URP = 10, LRP = 0, directory = os.getcwd(), inp_name = 'Input.inp', measure_time = 30):

    # URV = 5.01, LRV = 1.005 reverted to 5.248 and 1.054, respectively, 07JAN25
    with open(os.path.join(directory, inp_name), 'w') as f:
        print(f"\
*Pitot tube\n\
file={filename}.dat\n\
r/R={roverR}\n\
frequency=50000\n\
measuretime={measure_time}\n\
URV={URV}\n\
LRV={LRV}\n\
URP={URP}\n\
LRP={LRP}", file = f)
        
    return

    
def tdms_to_dat(infile:str, outfile = None):
    """_summary_

    :param infile: Input .tdms file to process
    :type infile: str
    :param outfile: Name of output .dat file, defaults to "infile".dat
    :type outfile: str, optional

    """
    from nptdms import TdmsFile

    tdms_file = TdmsFile.read(infile)
    data = []
    for group in tdms_file.groups():
        for channel in group.channels():
            data.append( channel[:] )

    if outfile is None:
        outfile = infile.strip('tdms') + 'dat'
    
    data = np.asarray(data)

    np.savetxt(outfile, data.T)

    return

def process_dir(target_dir:str, probe_number:str, r01:float, r02:float, r03:float, r12:float, r13:float, r23:float, roverR = None, measure_time = 30,
                signal_output=0, detailed_output=0, multiprocess = False, num_cpus = None, mode = "probe"):
    """ Runs MIDAS for every dat file in a given directory

    Makes a new folder, auto_reprocessed_data_TIMESTAMP, where the .tab files will be put.

    Inputs:
     - target_dir, directory to process
     - Probe number, for identification
     - Probe measurements (r01, r02, etc.). In mm. Same as in .inp file. For Pitot tube, these are repurposed as URV, LRV, URP, LRP
     - signal_output, makes the _MedianSig, _NormSig, etc. files
     - mode, either "probe" or "pitot"

    Outputs:
     - Returns name of directory the reprocessed files are in
    
    """
    os.chdir(target_dir)
    current_time = datetime.now()
    timestamp = f"{current_time.month}-{current_time.day}-{current_time.year}_{current_time.hour}-{current_time.minute}"

    reprocessed_dir = os.path.join(target_dir, 'auto_reprocessed_data_'+ timestamp)

    if not os.path.isdir( reprocessed_dir ):
        os.makedirs( reprocessed_dir )

    if mode == 'probe':
        try:
            copy2(os.path.join(target_dir, "MIDASv1.14d.exe"), reprocessed_dir)
        except FileNotFoundError:
            copy2(os.path.join("Z:\TRSL\PITA", "MIDASv1.14d.exe"), reprocessed_dir)
    elif mode == 'pitot':
        try:
            copy2(os.path.join(target_dir, "PPv1.exe"), reprocessed_dir)
        except FileNotFoundError:
            copy2(os.path.join("Z:\TRSL\PITA", "PPv1.exe"), reprocessed_dir)
    

    os.chdir(reprocessed_dir)

    if not multiprocess:
        for file in os.listdir(target_dir):
            # print(file)
            if file.split('.')[-1] == 'dat':
                copy2(os.path.join(target_dir, file), reprocessed_dir)

                if roverR is None:
                    roverR = 0.1 * int(file[1])
                
                if mode == 'probe':
                    write_inp(roverR, file.replace('.dat', ''), probe_number = probe_number, 
                              r01=r01, r02=r02, r03=r03, r12=r12, r13=r13, r23=r23, 
                          directory=reprocessed_dir, signalOutput=signal_output, detailedOutput=detailed_output, measure_time=measure_time
                          )
                    comp_process = run(os.path.join(reprocessed_dir, 'MIDASv1.14d.exe'), cwd = reprocessed_dir, shell=True)
                elif mode == 'pitot':
                    write_pitot_inp(roverR, file.replace('.dat', ''), measure_time=measure_time,
                                    URV = r01, LRV = r02,
                                    URP= r03, LRP = r12,
                                    directory=reprocessed_dir
                                    )
                    comp_process = run(os.path.join(reprocessed_dir, 'PPv1.exe'), cwd = reprocessed_dir, shell=True)

                if comp_process.returncode != 0:
                    print(comp_process)
                try:
                    os.remove(os.path.join(reprocessed_dir, file))
                except OSError as e:
                    print("Failed to remove .dat file, ", e)

    else:
        import multiprocessing
        import subprocess
        from multiprocessing.pool import ThreadPool
        
        def write_inp_run_midas(file, roverR):

            input_name = file.strip('.dat') + "_input.inp"

            if roverR is None:
                roverR = 0.1 * int(file[1])
            
            if mode == 'probe':
                write_inp(roverR, file.replace('.dat', ''), probe_number = probe_number, 
                          r01=r01, r02=r02, r03=r03, r12=r12, r13=r13, r23=r23, 
                          directory=reprocessed_dir, signalOutput=signal_output, detailedOutput=detailed_output, inp_name=input_name, measure_time=measure_time)
                p = subprocess.Popen(["MIDASv1.14d.exe", input_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            elif mode == 'pitot':
                write_pitot_inp(roverR, file.replace('.dat', ''), 
                                URV = r01, LRV = r02,
                                URP= r03, LRP = r12,
                                directory=reprocessed_dir,
                                inp_name=input_name, measure_time=measure_time)
                p = subprocess.Popen(["PPv1.exe", input_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            out, err = p.communicate()
            return (out, err)
        
        if num_cpus is None:
            num_cpus = multiprocessing.cpu_count()
        pool = ThreadPool(num_cpus)
        results = []

        for file in os.listdir(target_dir):
            if file.split('.')[-1] == 'dat':
                # print(file)
                if (mode == 'probe') and ('pitot' not in file):
                    copy2(os.path.join(target_dir, file), reprocessed_dir)
                    
                    results.append(pool.apply_async(write_inp_run_midas, (file, roverR) ) )
                
                elif (mode == 'pitot') and ('pitot' in file):
                    copy2(os.path.join(target_dir, file), reprocessed_dir)
                    
                    results.append(pool.apply_async(write_inp_run_midas, (file, roverR) ) )

        pool.close()
        pool.join()

        for result in results:
            out, err = result.get()
            print("out: {} err: {}".format(out, err))

        # for file in os.listdir(target_dir):
        #     try:
        #         os.remove(os.path.join(reprocessed_dir, file))
        #     except OSError as e:
        #         print("Failed to remove .dat file, ", e)

    return reprocessed_dir