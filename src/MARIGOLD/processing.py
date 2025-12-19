from .config import *

"""
This file contains functions related to interacting with raw data files from the NI card, MIDAS or PP

- write_inp, for MIDAS
- write_pitot_inp, for PP
- tdms_to_dat, converts an NI .tdms file to a .dat file
- process_dir, processes entire directory

"""

def write_inp(roverR, filename, probe_number = 'Dummy', probe_type = 4, set_vals = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], directory = os.getcwd(), detailedOutput=0, signalOutput=0, inp_name = 'Input.inp', measure_time = 30):
    """ Write an .inp file for MIDAS
    
    """
    # Unpack probe dimensions
    [r01, r02, r03, r12, r13, r23] = set_vals

    with open(os.path.join(directory, inp_name), 'w') as f:
        f.write(

f"*PROBE NUMBER: {probe_number}\n\
probetype={probe_type}\n\
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

def write_pitot_inp(roverR, filename, set_vals, directory = os.getcwd(), inp_name = 'Input.inp', measure_time = 30):

    # Unpack pressure bounds
    [URV, LRV, URP, LRP, _, _] = set_vals

    with open(os.path.join(directory, inp_name), 'w') as f:
        print(

f"*Pitot tube\n\
file={filename}.dat\n\
r/R={roverR}\n\
frequency=50000\n\
measuretime={measure_time}\n\
URV={URV}\n\
LRV={LRV}\n\
URP={URP}\n\
LRP={LRP}", file = f

)
        
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

def process_dir(target_dir:str, probe_number:str, set_vals,
                probe_type = 4, roverR = None, measure_time = 30,
                signal_output = 0, detailed_output = 0,
                multiprocess = False, num_cpus = None,
                mode = "probe", midas = "MIDASv1.14d.exe"):
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

    # Change directory to target
    os.chdir(target_dir)

    # Make reprocessed directory
    current_time = datetime.now()
    timestamp = f"{current_time.month}-{current_time.day}-{current_time.year}_{current_time.hour}-{current_time.minute}"

    reprocessed_dir = os.path.join(target_dir, 'auto_reprocessed_data_'+ timestamp)

    if not os.path.isdir( reprocessed_dir ):
        os.makedirs( reprocessed_dir )

    # Copy executable to reprocessed directory
    if mode == 'probe':
        copy2(os.path.join(target_dir, midas), reprocessed_dir)

    elif mode == 'pitot':
        copy2(os.path.join(target_dir, "PPv1.exe"), reprocessed_dir)    

    os.chdir(reprocessed_dir)

    # Loop through files
    if not multiprocess:
        for file in os.listdir(target_dir):
            if file.split('.')[-1] == 'dat':
                # Copy query file
                copy2(os.path.join(target_dir, file), reprocessed_dir)

                # Auto-detect radial location
                if roverR is None:
                    roverR = 0.1 * int(file[1])
                
                # Write input and run executable
                if mode == 'probe':
                    write_inp(roverR, file.replace('.dat', ''), probe_number = probe_number, probe_type = probe_type, set_vals = set_vals, 
                        directory = reprocessed_dir, signalOutput = signal_output, detailedOutput = detailed_output, measure_time=measure_time
                        )
                    
                    comp_process = run(os.path.join(reprocessed_dir, midas), cwd = reprocessed_dir, shell=True)
                    
                elif mode == 'pitot':
                    write_pitot_inp(roverR, file.replace('.dat', ''), measure_time = measure_time, set_vals = set_vals,
                        directory = reprocessed_dir
                        )
                    
                    comp_process = run(os.path.join(reprocessed_dir, 'PPv1.exe'), cwd = reprocessed_dir, shell=True)

                if comp_process.returncode != 0:
                    print(comp_process)
                
                # Delete copied file
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
                write_inp(roverR, file.replace('.dat', ''), probe_number = probe_number, probe_type = probe_type, set_vals = set_vals, 
                    directory = reprocessed_dir, signalOutput = signal_output, detailedOutput = detailed_output, inp_name = input_name, measure_time = measure_time
                    )
                
                p = subprocess.Popen([midas, input_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            elif mode == 'pitot':
                write_pitot_inp(roverR, file.replace('.dat', ''), set_vals = set_vals,
                    directory = reprocessed_dir, inp_name = input_name, measure_time = measure_time
                    )
                
                p = subprocess.Popen(["PPv1.exe", input_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            out, err = p.communicate()
            return (out, err)
        
        if num_cpus is None:
            num_cpus = multiprocessing.cpu_count()

        pool = ThreadPool(num_cpus)
        results = []

        for file in os.listdir(target_dir):
            if file.split('.')[-1] == 'dat':
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

    return reprocessed_dir