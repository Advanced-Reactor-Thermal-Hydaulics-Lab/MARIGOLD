
from .config import *
from .Condition import Condition
from .Iskandrani_Condition import Iskandrani_Condition
from .Yang_Condition import Yang_Condition

import re
import xlrd

def extractProbeData(dump_file = 'database.dat', in_dir = [], require_terms = None, skip_terms = ['CFD', 'Copy'],
                     extract_Ryan = True, Ryan_path = 'Z:\\TRSL\\PITA\\Data\\LocalData\\spreadsheets\\PITA',
                     extract_Kong = True, Kong_path = 'Z:\\TRSL\\PITA\\Data\\LocalData\\spreadsheets\\101.6mm',
                     extract_Talley = True, Talley_path = 'Z:\\TRSL\\PITA\\Data\\LocalData\\spreadsheets\\38.1mm',
                     extract_Yadav = True, Yadav_path = 'Z:\\TRSL\\PITA\\Data\\LocalData\\spreadsheets\\50.8mm'
                     ) -> None:
    debug = True
    all_conditions = []

    ### PITA Data ###

    if extract_Ryan:
        path = Ryan_path

        Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 33)], [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)], [i for i in range(208, 233)] ]))
        Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)] ]))

        #print(Q1_ranges, Q2_ranges)

        for file in os.listdir(path):
            if debug: print(file, file=debugFID)
            if file.split('.')[-1] == 'xlsx':
                
                # Check if the file has any skipped/required terms
                if any(term in file for term in skip_terms):
                    if debug: print(f"Skipping {file}", file=debugFID)
                    continue

                if require_terms:
                    if all(term not in file for term in require_terms):
                        if debug: print(f"Skipping {file}", file=debugFID)
                        continue
                
                #if debug: print(path, file=debugFID)
                
                try:
                    wb = op.load_workbook(filename=os.path.join(path, file), data_only=True)
                except:
                    print(f"Error reading wb: {file}\nSkipping...")
                    continue
                
                try:
                    jf = float(file.split('_')[1].strip('jf'))
                    jgref = float(file.split('_')[2].strip('jg'))
                    port = file.split('_')[3].strip('.xlsx')
                    theta = float(file.split('_')[0].strip('deg'))
                except:
                    print(f'Warning: Non-standard excel file name {file}. Skipping...')
                    continue

                ws = wb['1']
                jgloc = ws['U23'].value

                newCond = Condition(jgref, jgloc, jf, theta, port, 'Ryan')

                if newCond not in all_conditions:
                    all_conditions.append(newCond)
                    cond = newCond
                else:
                    cond = all_conditions[ all_conditions.index(newCond) ]
                
                ws = wb['2']

                cond.area_avg_void_sheet = ws['G266'].value
                
                for phi, indices in Q1_ranges:
                    first_phi = True
                    for i in indices:
                        if ws[f'K{i}'].value:

                            try:
                                roverR = float(ws[f'A{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...', file=debugFID)
                                continue
                            
                            midas_output = []
                            for cell in ws[f'A{i}':f'BD{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output", file=debugFID)
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)

                            try:
                                cond.data[phi].update({roverR: data})
                            except KeyError as e:
                                if first_phi:
                                    pass
                                else:
                                    if debug: print("Not my first phi, for some reaseon", e, file=debugFID)
                                cond.data.update( {phi:{}} )
                                first_phi = False
                                cond.data[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.data[phi].update({roverR: data})


                for phi, indices in Q2_ranges:
                    for i in indices:
                        if ws[f'DC{i}'].value:

                            try:
                                roverR = float(ws[f'CO{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...')
                                continue
                            
                            midas_output = []
                            for cell in ws[f'CO{i}':f'ER{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                print("Warning, tab_keys not the same length as midas_output")
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)

                            try:
                                cond.data[phi].update({roverR: data})
                            except KeyError:
                                cond.data.update( {phi:{}} )
                                cond.data[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.data[phi].update({roverR: data})

    
    ### 38.1 mm data (Talley) ###

    if extract_Talley:
        path = Talley_path
        #print(Q1_ranges, Q2_ranges)
        phis = [90, 67.5, 45, 22.5, 0]

        for file in os.listdir(path):
            if debug: print(file, file=debugFID)
            if file.split('.')[-1] == 'xlsx':
                
                # Check if the file has any skipped/required terms
                if any(term in file for term in skip_terms):
                    if debug: print(f"Skipping {file}", file=debugFID)
                    continue

                if any(term not in file for term in require_terms):
                    if debug: print(f"Skipping {file}", file=debugFID)
                    continue
                
                #if debug: print(path, file=debugFID)
                
                try:
                    wb = op.load_workbook(filename=os.path.join(path, file), data_only=True)
                except:
                    print(f"Error reading wb: {file}\nSkipping...")
                    continue
                
                try:
                    jf = float(file.split('_')[1].strip('jf'))
                    jgref = float(file.split('_')[2].strip('jg'))
                    port = file.split('_')[3].strip('.xlsx')
                    theta = float(file.split('_')[0].strip('deg'))
                except:
                    print(f'Warning: Non-standard excel file name {file}. Skipping...')
                    continue
                
                try:
                    ws = wb['2']
                    jgloc = ws['C2'].value
                    old = False
                except:
                    #print(f"Warning: Old format file {file}")
                    ws = wb['Sheet1']
                    jgloc = ws['C3'].value
                    old = True
                

                newCond = Condition(jgref, jgloc, jf, theta, port, 'Talley')

                if newCond not in all_conditions:
                    all_conditions.append(newCond)
                    cond = newCond
                else:
                    cond = all_conditions[ all_conditions.index(newCond) ]
                
                i = 0
                
                phi_counter = 0
                next = False
                while phi_counter < 5:
                    i += 1

                    if ws[f'E{i}'].value == 'Spherical' or  ws[f'C{i}'].value == 'Spherical':
                        if debug: print(f'found header in row {i}', file=debugFID)
                        next = True
                        continue

                    if next:
                        try:
                            # Hopefully in the part of the sheet with data
                            roverR = float(ws[f'A{i}'].value)
                        except:
                            # Done reading data
                            phi_counter += 1
                            next = False
                            continue


                        midas_output = []
                        data = deepcopy(zero_data)

                        if old and ws[f'F{i}'].value:
                            # use old tab keys
                            for cell in ws[f'A{i}':f'AA{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(old_tab_keys) == len( midas_output ):
                                data = dict( zip( old_tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, old tab_keys not the same length as midas_output")
                                data = dict( zip( old_tab_keys, midas_output ))
                                if debug:
                                    print("old tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)
                        elif ws[f'K{i}'].value:
                            # use new tab keys
                            roverR = float(ws[f'A{i}'].value)
                            midas_output = []
                            for cell in ws[f'A{i}':f'BD{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output")
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)

                        phi = phis[phi_counter]
                        try:
                            cond.data[phi].update({roverR: data})
                        except KeyError:
                            cond.data.update( {phi:{}} )
                            cond.data[phi].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.data[phi].update({roverR: data})

                #print('total rows read:', i)

    ### 101.6mm Data ###

    if extract_Kong:
        path = Kong_path
        Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 33)], [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)], [i for i in range(208, 233)] ]))
        Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)] ]))

        #print(Q1_ranges, Q2_ranges)

        for file in os.listdir(path):
            if debug: print(file, file=debugFID)
            if file.split('.')[-1] == 'xlsx':
                
                # Check if the file has any skipped/required terms
                if any(term in file for term in skip_terms):
                    if debug: print(f"Skipping {file}", file=debugFID)
                    continue

                if any(term not in file for term in require_terms):
                    if debug: print(f"Skipping {file}", file=debugFID)
                    continue
                
                #if debug: print(path, file=debugFID)
                
                try:
                    wb = op.load_workbook(filename=os.path.join(path, file), data_only=True)
                except:
                    print(f"Error reading wb: {file}\nSkipping...")
                    continue
                
                try:
                    jf = float(file.split('_')[1].strip('jf'))
                    jgref = float(file.split('_')[2].strip('jg'))
                    port = file.split('_')[3].strip('.xlsx')
                    theta = float(file.split('_')[0].strip('deg'))
                except:
                    print(f'Warning: Non-standard excel file name {file}. Skipping...')
                    continue
                
                try:
                    ws = wb['1']
                    jgloc = ws['U23'].value
                except:
                    if debug: print(f'Warning: no sheet "1" found for {file}. Using sheet "2" C2 for jgloc', file=debugFID)
                    ws = wb['2']
                    jgloc = ws['C2'].value
                

                newCond = Condition(jgref, jgloc, jf, theta, port, 'Kong')

                if newCond not in all_conditions:
                    all_conditions.append(newCond)
                    cond = newCond
                else:
                    cond = all_conditions[ all_conditions.index(newCond) ]
                
                ws = wb['2']
                cond.area_avg_void_sheet = ws['G266'].value
                
                for phi, indices in Q1_ranges:
                    for i in indices:
                        if ws[f'K{i}'].value:

                            try:
                                roverR = float(ws[f'A{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...')
                                continue
                            
                            midas_output = []
                            for cell in ws[f'A{i}':f'BD{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output")
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)

                            try:
                                cond.data[phi].update({roverR: data})
                            except KeyError:
                                cond.data.update( {phi:{}} )
                                cond.data[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.data[phi].update({roverR: data})


                for phi, indices in Q2_ranges:
                    for i in indices:
                        if ws[f'DC{i}'].value:

                            try:
                                roverR = float(ws[f'CN{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...')
                                continue
                            
                            midas_output = []
                            for cell in ws[f'CN{i}':f'EQ{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output")
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)

                            try:
                                cond.data[phi].update({roverR: data})
                            except KeyError:
                                cond.data.update( {phi:{}} )
                                cond.data[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.data[phi].update({roverR: data})
                
    
    ### Yadav Data ###

    if extract_Yadav:
        path = Yadav_path
        print("Reading Yadav's data")
        Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 29)], [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)], [i for i in range(192, 213)] ]))
        Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)] ]))

        for file in os.listdir(path):
            if debug: print(file, file=debugFID)
            if file.split('.')[-1] == 'xlsx':
                
                # Check if the file has any skipped/required terms
                if any(term in file for term in skip_terms):
                    if debug: print(f"Skipping {file}", file=debugFID)
                    continue

                if all(term not in file for term in require_terms):
                    if debug: print(f"Skipping {file}", file=debugFID)
                    continue
                
                #if debug: print(path, file=debugFID)
                
                try:
                    wb = op.load_workbook(filename=os.path.join(path, file), data_only=True)
                except:
                    print(f"Error reading wb: {file}\nSkipping...")
                    continue
                
                try:
                    jf = float(file.split('_')[1].strip('jf'))
                    jgref = float(file.split('_')[2].strip('jg'))
                    port = file.split('_')[3].strip('.xlsx')
                    theta = float(file.split('_')[0].strip('deg'))
                except:
                    print(f'Warning: Non-standard excel file name {file}. Skipping...')
                    continue

                ws = wb['1']

                # jglocs = ['O16', 'O17', 'O18', 'O19', 'O20']
                # jgloc = ws[jglocs[int(port.strip('P'))-1]].value
                jgloc = jgref
                
                newCond = Condition(jgref, jgloc, jf, theta, port, 'Yadav')

                if newCond not in all_conditions:
                    all_conditions.append(newCond)
                    cond = newCond
                else:
                    cond = all_conditions[ all_conditions.index(newCond) ]
                
                ws = wb['2']

                cond.area_avg_void_sheet = ws['G246'].value
                cond.jgloc = ws['C3'].value
                
                for phi, indices in Q1_ranges:
                    for i in indices:
                        if ws[f'AA{i}'].value:
                            if (ws[f'A{i}'].value == 0) and (ws[f'AA{i}'].value == 0.0064) and (ws[f'AD{i}'].value == 13.45) and (ws[f'AY{i}'].value == 1.682):
                                # It's dummy data, skip it
                                if debug: print(f'Dummy data found in row {i} in sheet {file}. Skipping...', file=debugFID)
                                continue
                            
                            try:
                                roverR = float(ws[f'A{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...', file=debugFID)
                                continue
                            
                            midas_output = []
                            for cell in ws[f'A{i}':f'BD{i}'][0]:
                                if cell.value == '---':
                                    midas_output.append(0)
                                elif cell.value is None:
                                    if debug: print(f"None found in {cell}")
                                    midas_output.append(0)
                                else:
                                    midas_output.append(float(cell.value))
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output", file=debugFID)
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)

                            try:
                                cond.data[phi].update({roverR: data})
                            except KeyError:
                                cond.data.update( {phi:{}} )
                                cond.data[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.data[phi].update({roverR: data})


            for phi, indices in Q2_ranges:
                for i in indices:
                    if ws[f'DR{i}'].value:

                        try:
                            roverR = float(ws[f'CP{i}'].value)
                        except:
                            if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...')
                            continue
                        
                        midas_output = []
                        for cell in ws[f'CQ{i}':f'ES{i}'][0]:
                            if cell.value == '---':
                                midas_output.append(0)
                            elif cell.value is None:
                                if debug: print(f"None found in {cell}")
                                midas_output.append(0)
                            else:
                                midas_output.append(float(cell.value))
                        
                        if len(tab_keys) == len( midas_output ):
                            data = dict( zip( tab_keys, midas_output ))
                        else:
                            if debug: print("Warning, tab_keys not the same length as midas_output")
                            data = dict( zip( tab_keys, midas_output ))
                            if debug:
                                print("tab_keys not the same length as midas_output", file=debugFID)
                                print(tab_keys, midas_output, file=debugFID)

                        try:
                            cond.data[phi].update({roverR: data})
                        except KeyError:
                            cond.data.update( {phi:{}} )
                            cond.data[phi].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.data[phi].update({roverR: data}) 

    if debug and False:
        for cond in all_conditions:
            cond.pretty_print()
    
    with open(dump_file, 'wb') as g:
        pickle.dump(all_conditions, g)
    return


def extractLocalDataFromDir(path:str, dump_file = 'database.dat', in_dir = [], require_terms = ['jf'], 
                            skip_terms = ['CFD', 'Copy'], sheet_type = 'adix_template', append_to_json = None,
                            pitot_sheet = False, print_sheets = False,
                            **kwargs) -> None:
    
    """Function for getting all local data from spreadsheets in a directory, path

       Does not recursively descend, only checks in the path given
    
       Still under construction, but should support sheet types
        - 'adix_template4'
        - 'ryan_template'
        - 'adix_template' (maybe rename this quan_template)
        - 'bettis_template'
        - 'talley_template'
        - 'yadav_template'

       Also can try to infer the sheet type, xlsm will be adix_template4, if it has P5, 6, or 7 it will
       be classified as an adix_template, different angles ryan, etc. The inference can also be made by 
       appending _adix or _quan or _ryan to the excel sheets being processed

       Custom sheet types may be supported, but they still have to generally follow the classic template
       structure. The setup information has to be in sheet '1' and all the local data in sheet '2'. Ranges
       must be specified as a list of lists

       Q1_ranges = [[angle1, [index1, index2.,..]], [angle2, [index1, index2.,..]], ...]

       with the starts and ends

       [Q1_start, Q1_end, Q2_start, ...]

       Add on pitot starts/ends if pitot_sheet = True

       Need to specify Q1_check, Q2_check, etc. as well

    """
    all_conditions = []

    for file in os.listdir(path):
    
        # print(file)
        if print_sheets: print(file)

        if file.split('.')[-1] == 'xlsx' or file.split('.')[-1] == 'xlsm' or file.split('.')[-1] == 'xls':

            # Check if the file has any skipped/required terms
            if any(term in file for term in skip_terms):
                if debug: print(f"Skipping {file}", file=debugFID)
                continue

            # if all(term not in file for term in require_terms):
            #     if debug: print(f"Skipping {file}", file=debugFID)
            #     continue
            
            # if debug: print(path, file=debugFID)
            
            try:
                if file.split('.')[-1] == 'xls':
                    wb = xlrd.open_workbook(filename=os.path.join(path, file))
                else:
                    wb = op.load_workbook(filename=os.path.join(path, file), data_only=True)                    
            except:
                print(f"Error reading wb: {file}\nSkipping...")
                continue
            
            try:
                jf = float(file.split('_')[1].strip('jf'))
                jgref = float(file.split('_')[2].strip('jg'))
                port = file.split('_')[3].strip('.xlsx').strip('.xlsm')
                theta = float(file.split('_')[0].strip('deg'))
            except:
                print(f'Warning: Non-standard excel file name {file}. Is this Bettis template?')
                
                if sheet_type == 'bettis_template' and 'Run' in file.split('_')[0]:
                    print("Yes, it is. Proceeding...")

                    theta = 90
                    port = file.split('_')[-1]
                    port_idx = int(''.join(re.findall(r'\d+',port)))
                
                    if file.split('_')[0] == 'Run1':
                        jf = 0.32
                        jgref = 0.047
                        
                        run_idx = 1
                    elif file.split('_')[0] == 'Run2':
                        jf = 0.95
                        jgref = 0.047

                        run_idx = 2
                    elif file.split('_')[0] == 'Run3':
                        jf = 1.89
                        jgref = 0.095
                        
                        run_idx = 3
                    elif file.split('_')[0] == 'Run4':
                        jf = 0.95
                        jgref = 0.187
                        
                        run_idx = 4
                    elif file.split('_')[0] == 'Run5':
                        jf = 1.89
                        jgref = 0.193
                        
                        run_idx = 5
                    elif file.split('_')[0] == 'Run6' and file.split('_')[1] == 'short':
                        jf = 0.63
                        jgref = 0.279
                        
                        run_idx = 6
                        
                        print("Ommitting short runs...")
                        continue
                    elif file.split('_')[0] == 'Run7' and file.split('_')[1] == 'short':
                        jf = 2.84
                        jgref = 0.287
                        
                        run_idx = 7

                        print("Ommitting short runs...")
                        continue
                    elif file.split('_')[0] == 'Run6' and file.split('_')[1] == 'long':
                        jf = 0.63
                        jgref = 0.279
                        
                        run_idx = 8
                    elif file.split('_')[0] == 'Run7' and file.split('_')[1] == 'long':
                        jf = 2.84
                        jgref = 0.287
                        
                        run_idx = 9
                    elif file.split('_')[0] == 'Run8':
                        jf = 1.89
                        jgref = 0.385
                        
                        run_idx = 10
                    elif file.split('_')[0] == 'Run9':
                        jf = 4.40
                        jgref = 0.940
                        
                        run_idx = 11
                    else:
                        print("Warning: Run number exceeds highest known run. Skipping...")
                        continue
                    
                else:
                    print(f'Nope. Skipping...')
                    continue
            # print(jf, jgref, port, theta)
            ############################################################################################################################
            #                                                                                                                          #
            #                                                       BETTIS DATA                                                        #
            #                                                                                                                          #
            ############################################################################################################################
            if sheet_type == 'bettis_template':

                '''
                # Sourced from Kim_Research_BubbleDoc\Bettis\OneGroupEvaluation\AllConditions
                jgloc_mat = [
                    [0.047 , 0.040 , 0.040 , 0.040 , None  , None  ],   # Run1
                    [0.047 , 0.040 , 0.040 , 0.040 , None  , None  ],   # Run2
                    [0.095 , 0.070 , 0.070 , 0.070 , None  , None  ],   # Run3
                    [0.187 , 0.140 , 0.150 , 0.150 , None  , None  ],   # Run4
                    [0.193 , 0.140 , 0.150 , 0.150 , None  , None  ],   # Run5
                    [0.279 , 0.220 , None  , 0.220 , None  , None  ],   # Run6_short
                    [0.287 , 0.210 , 0.220 , 0.220 , None  , None  ],   # Run7_short
                    [None  , 0.219 , None  , 0.234 , None  , 0.251 ],   # Run6_long
                    [None  , 0.212 , None  , None  , None  , 0.264 ],   # Run7_long
                    [None  , 0.288 , None  , 0.314 , None  , 0.346 ],   # Run8
                    [None  , 0.618 , None  , 0.716 , None  , 0.850 ],   # Run9
                ]

                alpha_mat = [
                    [0.066 , 0.089 , 0.080 , 0.083 , None  , None  ],   # Run1
                    [0.039 , 0.038 , 0.037 , 0.034 , None  , None  ],   # Run2
                    [0.035 , 0.038 , 0.031 , 0.034 , None  , None  ],   # Run3
                    [0.104 , 0.110 , 0.116 , 0.116 , None  , None  ],   # Run4
                    [0.064 , 0.059 , 0.062 , 0.059 , None  , None  ],   # Run5
                    [0.220 , 0.173 , None  , 0.248 , None  , None  ],   # Run6_short
                    [0.031 , 0.051 , None  , 0.058 , None  , None  ],   # Run7_short
                    [None  , 0.207 , None  , 0.222 , None  , 0.237 ],   # Run6_long
                    [None  , 0.067 , None  , None  , None  , 0.075 ],   # Run7_long
                    [None  , 0.121 , None  , 0.103 , None  , 0.124 ],   # Run8
                    [None  , 0.078 , None  , 0.095 , None  , 0.125 ]    # Run9
                ]

                ai_mat = [
                    [None  , 190.33, 175.51, 173.26, None  , None  ],   # Run1
                    [None  , 78.90 , 77.30 , 71.10 , None  , None  ],   # Run2
                    [None  , 87.80 , 76.70 , 78.60 , None  , None  ],   # Run3
                    [None  , 210.56, 223.89, 211.69, None  , None  ],   # Run4
                    [None  , 125.06, 135.80, 120.97, None  , None  ],   # Run5
                    [None  , 432.78, None  , 372.21, None  , None  ],   # Run6_short
                    [None  , 151.23, None  , 157.03, None  , None  ],   # Run7_short
                    [None  , 470.61, None  , 369.69, None  , 418.74],   # Run6_long
                    [None  , 199.34, None  , None  , None  , 222.64],   # Run7_long
                    [None  , 209.62, None  , 192.28, None  , 270.00],   # Run8
                    [None  , 188.53, None  , 258.80, None  , 345.16],   # Run9
                ]

                Dsm_mat = [
                    [None  , 2.56  , None  , 2.61  , None  , None  ],   # Run1
                    [None  , 2.87  , None  , 2.71  , None  , None  ],   # Run2
                    [None  , 2.55  , None  , 2.59  , None  , None  ],   # Run3
                    [None  , 3.06  , None  , 3.06  , None  , None  ],   # Run4
                    [None  , 2.92  , None  , 2.92  , None  , None  ],   # Run5
                    [None  , 2.33  , None  , 3.81  , None  , None  ],   # Run6_short
                    [None  , 2.03  , None  , 2.21  , None  , None  ],   # Run7_short
                    [None  , 2.61  , None  , 3.54  , None  , 3.28  ],   # Run6_long
                    [None  , 1.98  , None  , None  , None  , 2.01  ],   # Run7_long, last one at 140.1 L/D
                    [None  , 2.67  , None  , 3.26  , None  , 2.78  ],   # Run8
                    [None  , 2.43  , None  , 2.19  , None  , 2.18  ],   # Run9
                ]

                # Sourced from thesis Table A.E.(a), and MATLAB initcond.m
                pz_mat = [
                    [None  , 30603.60 , 27710.01 , 24816.42 , None  , None    ],  # Run1
                    [None  , 32735.69 , 29567.55 , 26399.42 , None  , None    ],  # Run2
                    [None  , 35073.44 , 31433.39 , 27793.33 , None  , None    ],  # Run3
                    [None  , 31486.61 , 28501.50 , 25516.40 , None  , None    ],  # Run4
                    [None  , 38464.37 , 34860.21 , 31256.04 , None  , None    ],  # Run5
                    [None  , 29584.83 , 27028.58 , 24472.33 , None  , None    ],  # Run6_short
                    [None  , 36275.59 , 32140.29 , 28004.99 , None  , None    ],  # Run7_short
                    [None  , 28213.762, None     , 19814.36 , None  , 11415   ],  # Run6_long
                    [None  , 40150.62 , None     , 25970.05 , None  , 11789.5 ],  # Run7_long
                    [None  , 34289.443, None     , 22869.15 , None  , 11448.8 ],  # Run8
                    [None  , 52254.526, None     , 31247.62 , None  , 10240.7 ],  # Run9
                ]

                LoverD_mat = [8.02, 34.76, 61.49, 88.22, 114.96, 141.70]
                dpdz_mat = [9250.38, 10128.04, 11636.70, 9542.94, 11521.98, 8171.95, 13219.94, 8267.13, 13957.25, 11240.45, 20676.09]
                '''

                # Sourced from \Kim_ThesisFolder\Thesis_Original\OneGroupEvaluation\AllConditions
                jgloc_mat = [
                    [0.047 , 0.040 , 0.040 , 0.040 , None  , None  ],   # Run1
                    [0.047 , 0.040 , 0.040 , 0.040 , None  , None  ],   # Run2
                    [0.095 , 0.070 , 0.070 , 0.070 , None  , None  ],   # Run3
                    [0.187 , 0.140 , 0.150 , 0.150 , None  , None  ],   # Run4
                    [0.193 , 0.140 , 0.150 , 0.150 , None  , None  ],   # Run5
                    [0.279 , 0.220 , None  , 0.220 , None  , None  ],   # Run6_short
                    [0.287 , 0.210 , 0.220 , 0.220 , None  , None  ],   # Run7_short
                    [None  , 0.28  , None  , 0.234 , None  , 0.251 ],   # Run6_long
                    [None  , 0.3   , None  , None  , None  , 0.264 ],   # Run7_long
                    [None  , 0.385 , None  , 0.314 , None  , 0.346 ],   # Run8
                    [None  , 0.94  , None  , 0.716 , None  , 0.850 ],   # Run9
                ]

                alpha_mat = [
                    [0.066 , 0.0839, 0.080 , 0.083 , None  , None  ],   # Run1
                    [0.039 , 0.038 , 0.037 , 0.034 , None  , None  ],   # Run2
                    [0.035 , 0.038 , 0.031 , 0.034 , None  , None  ],   # Run3
                    [0.104 , 0.1097, 0.116 , 0.116 , None  , None  ],   # Run4
                    [0.064 , 0.0592, 0.062 , 0.059 , None  , None  ],   # Run5
                    [0.220 , 0.1727, None  , 0.248 , None  , None  ],   # Run6_short
                    [0.031 , 0.051 , None  , 0.058 , None  , None  ],   # Run7_short
                    [None  , 0.2068, None  , 0.222 , None  , 0.237 ],   # Run6_long
                    [None  , 0.0674, None  , None  , None  , 0.075 ],   # Run7_long
                    [None  , 0.1149, None  , 0.103 , None  , 0.124 ],   # Run8
                    [None  , 0.0782, None  , 0.095 , None  , 0.125 ]    # Run9
                ]

                ai_mat = [
                    [None  , 190.33, 175.51, 173.26, None  , None  ],   # Run1
                    [None  , 82.28 , 77.30 , 71.10 , None  , None  ],   # Run2
                    [None  , 87.70 , 76.70 , 78.60 , None  , None  ],   # Run3
                    [None  , 210.56, 223.89, 211.69, None  , None  ],   # Run4
                    [None  , 125.06, 135.80, 120.97, None  , None  ],   # Run5
                    [None  , 432.78, None  , 372.21, None  , None  ],   # Run6_short
                    [None  , 151.23, None  , 157.03, None  , None  ],   # Run7_short
                    [None  , 470.61, None  , 369.69, None  , 418.74],   # Run6_long
                    [None  , 199.34, None  , None  , None  , 222.64],   # Run7_long
                    [None  , 209.62, None  , 192.28, None  , 270.00],   # Run8
                    [None  , 188.53, None  , 258.80, None  , 345.16],   # Run9
                ]

                Dsm_mat = [
                    [None  , 2.56  , None  , 2.61  , None  , None  ],   # Run1
                    [None  , 2.87  , None  , 2.71  , None  , None  ],   # Run2
                    [None  , 2.55  , None  , 2.59  , None  , None  ],   # Run3
                    [None  , 3.06  , None  , 3.06  , None  , None  ],   # Run4
                    [None  , 2.92  , None  , 2.92  , None  , None  ],   # Run5
                    [None  , 2.33  , None  , 3.81  , None  , None  ],   # Run6_short
                    [None  , 2.03  , None  , 2.21  , None  , None  ],   # Run7_short
                    [None  , 2.61  , None  , 3.54  , None  , 3.28  ],   # Run6_long
                    [None  , 1.98  , None  , None  , None  , 2.01  ],   # Run7_long, last one at 140.1 L/D
                    [None  , 2.67  , None  , 3.26  , None  , 2.78  ],   # Run8
                    [None  , 2.43  , None  , 2.19  , None  , 2.18  ],   # Run9
                ]

                # Sourced from thesis Table A.E.(a), and MATLAB initcond.m
                pz_mat = [
                    [None  , 30603.60 , 27710.01 , 24816.42 , None  , None    ],  # Run1
                    [None  , 32735.69 , 29567.55 , 26399.42 , None  , None    ],  # Run2
                    [None  , 35073.44 , 31433.39 , 27793.33 , None  , None    ],  # Run3
                    [None  , 31486.61 , 28501.50 , 25516.40 , None  , None    ],  # Run4
                    [None  , 38464.37 , 34860.21 , 31256.04 , None  , None    ],  # Run5
                    [None  , 29584.83 , 27028.58 , 24472.33 , None  , None    ],  # Run6_short
                    [None  , 36275.59 , 32140.29 , 28004.99 , None  , None    ],  # Run7_short
                    [None  , 28213.762, None     , 19814.36 , None  , 11415   ],  # Run6_long
                    [None  , 40150.62 , None     , 25970.05 , None  , 11789.5 ],  # Run7_long
                    [None  , 34289.443, None     , 22869.15 , None  , 11448.8 ],  # Run8
                    [None  , 52254.526, None     , 31247.62 , None  , 10240.7 ],  # Run9
                ]

                LoverD_mat = [8.02, 34.76, 61.49, 88.22, 114.96, 141.70]
                dpdz_mat = [9250.38, 10128.04, 11636.70, 9542.94, 11521.98, 8171.95, 13219.94, 8267.13, 13957.25, 11240.45, 20676.09]


                jgloc                   = jgloc_mat[run_idx-1][port_idx-1]
                area_avg_void_sheet     = alpha_mat[run_idx-1][port_idx-1]
                area_avg_ai_sheet       = ai_mat[run_idx-1][port_idx-1]
                area_avg_Dsm_sheet      = Dsm_mat[run_idx-1][port_idx-1]
                pz                      = pz_mat[run_idx-1][port_idx-1] + 101330
                LoverD                  = LoverD_mat[port_idx-1]
                dpdz                    = dpdz_mat[run_idx-1]
                
                newCond = Condition(jgref, jgloc, jf, theta, port, sheet_type.split('_')[0])

                if newCond not in all_conditions:
                    all_conditions.append(newCond)
                    cond = newCond
                else:
                    cond = all_conditions[ all_conditions.index(newCond) ]

                cond.Dh = 4 * 0.20 * 0.01 / 2 / (0.20 + 0.01)        # 1 x 20 cm^2 rectangular channel

                cond.area_avg_void_sheet = area_avg_void_sheet
                cond.area_avg_ai_sheet = area_avg_ai_sheet
                cond.area_avg_Dsm_sheet = area_avg_Dsm_sheet
                cond.pz = pz
                cond.LoverD = LoverD
                cond.dpdz = dpdz

                cond.jgatm = jgloc * pz / 101330        # It's important that P_atm is 101330 if you want to match with old results
            
            ############################################################################################################################
            #                                                                                                                          #
            #                                                       TALLEY DATA                                                        #
            #                                                                                                                          #
            ############################################################################################################################
            elif sheet_type == 'talley_template':
                #print(Q1_ranges, Q2_ranges)
                
                ws = wb['Sheet1']
                jgatm = ws['C3'].value
                old = True

                jgloc = jgref
                
                newCond = Condition(jgatm, jgloc, jf, theta, port, 'Talley')

                if newCond not in all_conditions:
                    all_conditions.append(newCond)
                    cond = newCond
                else:
                    cond = all_conditions[ all_conditions.index(newCond) ]
                
                cond.jgatm = jgatm
                cond.jgloc = jgloc

                cond.area_avg_void_sheet = ws['B206'].value
                cond.area_avg_ai_sheet = ws['C206'].value
                
                ########################################################################################################################
                # Q1

                phis = [90, 67.5, 45, 22.5, 0]

                i = 0
                phi_counter = 0
                next = False
                while phi_counter < 5:
                    i += 1

                    if ws[f'E{i}'].value == 'Spherical' or  ws[f'C{i}'].value == 'Spherical':
                        if debug: print(f'found header in row {i}', file=debugFID)
                        next = True
                        continue

                    if next:
                        try:
                            # Hopefully in the part of the sheet with data
                            roverR = float(ws[f'A{i}'].value)
                        except:
                            # Done reading data
                            phi_counter += 1
                            next = False
                            continue

                        midas_output = []
                        data = deepcopy(zero_data)

                        if old and ws[f'F{i}'].value:
                            for cell in ws[f'A{i}':f'AA{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(old_tab_keys) == len( midas_output ):
                                data = dict( zip( old_tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, old tab_keys not the same length as midas_output")
                                data = dict( zip( old_tab_keys, midas_output ))
                                if debug:
                                    print("old tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)
                        elif ws[f'K{i}'].value:
                            roverR = float(ws[f'A{i}'].value)
                            midas_output = []
                            for cell in ws[f'A{i}':f'BD{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output")
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)

                        phi = phis[phi_counter]
                        try:
                            cond.data[phi].update({roverR: data})
                        except KeyError:
                            cond.data.update( {phi:{}} )
                            cond.data[phi].update({1.0: zero_data})
                            cond.data[phi].update({roverR: data})

                ########################################################################################################################
                # Q2
                phis = [112.5, 135, 157.5]

                i = 0
                phi_counter = 0
                next = False
                while phi_counter < 3 and i < 200:
                    i += 1

                    if ws[f'AS{i}'].value == 'Spherical':
                        if debug: print(f'found header in row {i}', file=debugFID)
                        next = True
                        continue

                    if next:
                        try:
                            # Hopefully in the part of the sheet with data
                            roverR = float(ws[f'AQ{i}'].value)
                        except:
                            # Done reading data
                            phi_counter += 1
                            next = False
                            continue

                        midas_output = []
                        data = deepcopy(zero_data)

                        if old and ws[f'AV{i}'].value:
                            for cell in ws[f'AQ{i}':f'BQ{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(old_tab_keys) == len( midas_output ):
                                data = dict( zip( old_tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, old tab_keys not the same length as midas_output")
                                data = dict( zip( old_tab_keys, midas_output ))
                                if debug:
                                    print("old tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)
                        
                        phi = phis[phi_counter]
                        try:
                            cond.data[phi].update({roverR: data})
                        except KeyError:
                            cond.data.update( {phi:{}} )
                            cond.data[phi].update({1.0: zero_data})
                            cond.data[phi].update({roverR: data})

            ############################################################################################################################
            #                                                                                                                          #
            #                                                       YADAV DATA                                                         #
            #                                                                                                                          #
            ############################################################################################################################
            elif sheet_type == 'yadav_template':

                potent_ranges = [ [i for i in range(8, 22)], [i for i in range(48, 62)], [i for i in range(87, 101)], [i for i in range(128, 142)], [i for i in range(168, 182)], [i for i in range(209, 223)], [i for i in range(250, 264)], [i for i in range(290, 304)] ]
                
                # ID is going to be at A1, A41, A80, A121, A161, A202, A243, A283
                # jgatm C3
                # jf C4

                ws = wb['1']

                jgatm = ws['C3'].value
                jgref = jgatm

                jglocs = []

                jgloc = jgatm       # Nope. Either need jgloc or P_loc, but I have neither -- how did Yadav get jgloc for his IATE script?
                
                if type(jgloc) == str:
                    raise TypeError(f"\n\tin {file}\n\tRead in jgloc: {ws['C3'].value}\n\tInvalid type for jgloc")

                '''
                8 cond x 7 port?

                % exp_z = [0.38, 1.75, 3.12, 3.59, 4.96, 8.16, 12.58];

                titles=['Run1: j_{g}=0.039 [m/s] & j_{f}=0.018 [m/s]'; 'Run2: j_{g}=0.039 [m/s] & j_{f}=0.682 [m/s]'...
                    ;'Run3: j_{g}=0.136 [m/s] & j_{f}=0.682 [m/s]';'Run4: j_{g}=0.138 [m/s] & j_{f}=2.336 [m/s]';...
                    'Run5: j_{g}=0.538 [m/s] & j_{f}=5.100 [m/s]';'Run6: j_{g}=1.234 [m/s] & j_{f}=5.100 [m/s]';...
                    'Run7: j_{g}=0.039 [m/s] & j_{f}=0.064 [m/s]';'Run8: j_{g}=0.039 [m/s] & j_{f}=0.020 [m/s]'];

                jg_exp = [
                    0.083	0.090	0.099	0.101	0.102	0.104	0.106
                    0.136	0.148	0.163	0.166	0.168	0.171	0.175	
                    0.200	0.217	0.237	0.241	0.244	0.248	0.254	
                    0.072	0.078	0.085	0.087	0.089	0.091	0.095	
                    0.118	0.128	0.140	0.143	0.145	0.149	0.156	
                    0.177	0.191	0.207	0.212	0.216	0.222	0.232
                    0.094	0.101	0.108	0.110	0.113	0.123	0.139
                    0.141	0.148	0.159	0.162	0.167	0.182	0.206
                    ]
                
                p_exp = [
                    65707       51705.5     37799.6     34874.5     33378.4     31206.5     28138.4
                    68051       54277.2     40722.9     37763.4     36191.4     34100.3     30710.1
                    72326       58793.3     45604.4     42651.8     40976.4     38958.3     35329.6
                    94527       79201.8     64227.2     60578.2     57737.5     53657.9     46885.2
                    97147       81952.8     67116.1     63577.4     60495.5     56277.9     49353.5
                    100801      85731.1     71121.9     67438.5     64391.0     59910.8     52876.7
                    148848      132598.1    115749.7	111996.5	106153.9	89964.0     68377.4
                    153108      136853.6    120051.4	116366.2	110189.3	92801.6     69618.0
                    ];
                '''

                newCond = Condition(jgref, jgloc, jf, theta, port, sheet_type.split('_')[0])

                if newCond not in all_conditions:
                    all_conditions.append(newCond)
                    cond = newCond
                else:
                    cond = all_conditions[ all_conditions.index(newCond) ]

                cond.jgatm = jgatm

                cond.area_avg_void_sheet = ws['J328'].value
                cond.area_avg_ai_sheet = ws['K328'].value

                for phi, indices in Q1_ranges:
                    for i in indices:
                        if ws[f'{Q1_check}{i}'].value:

                            try:
                                roverR = float(ws[f'{Q1_start}{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column {Q1_start} could not be floatified. Skipping...', file=debugFID)
                                continue
                            
                            midas_output = []
                            for cell in ws[f'{Q1_start}{i}':f'{Q1_end}{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output", file=debugFID)
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)
                            
                            try:
                                cond.data[phi].update({roverR: data})
                            except KeyError:
                                cond.data.update( {phi:{}} )
                                cond.data[phi].update({1.0: deepcopy(zero_data)})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.data[phi].update({roverR: data})

            ############################################################################################################################
            #                                                                                                                          #
            #                                                        PITA DATA                                                         #
            #                                                                                                                          #
            ############################################################################################################################
            # General PITA template structure holds
            else:
                if sheet_type == 'infer':
                    if file.split('.')[-1] == 'xlsm' or 'adix' in file:
                        sheet_type = 'adix_template4'
                    elif '60deg' in file or '30deg' in file or '80deg' in file or 'Ryan' in file or 'ryan' in file:
                        sheet_type = 'ryan_template'
                    elif 'P4' in file or 'P5' in file or 'P6' in file or 'P7' in file or 'Quan' in file or 'quan' in file:
                        sheet_type = 'adix_template' # Might rename this Quan template
                    else:
                        sheet_type = 'adix_template'

                if sheet_type == 'ryan_template':
                    Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 33)], [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)], [i for i in range(208, 233)] ]))
                    Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)] ]))
                    Q2_start = 'CP'
                    Q2_end = 'ES'
                    Q1_start = 'A'
                    Q1_end = 'BD'

                    Q1_check = 'K'
                    Q2_check = 'DA'
                
                elif sheet_type == 'adix_template':
                    Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 29)], [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)], [i for i in range(192, 213)] ]))
                    Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)] ]))
                    Q2_start = 'CP'
                    Q2_end = 'ES'
                    Q1_start = 'A'
                    Q1_end = 'BD'

                    Q1_check = 'K'
                    Q2_check = 'DA'
                
                elif sheet_type == 'adix_template4':
                    pitot_sheet = True
                    Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 29)], [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)], [i for i in range(192, 213)] ]))
                    Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)] ]))
                    Q2_start = 'CR'
                    Q2_end = 'ET'
                    Q1_start = 'A'
                    Q1_end = 'BD'
                    Q1_pitot_start = 'CF'
                    Q1_pitot_end = 'CO'
                    Q2_pitot_start = 'FV'
                    Q2_pitot_end = 'GE'

                    Q1_check = 'K'
                    Q2_check = 'DA'
                    Q1_pitot_check = 'CJ'
                    Q2_pitot_check = 'FZ'

                elif sheet_type == 'custom' or sheet_type == 'Custom':
                    print('Hopefully you specified all the ranges, starts, ends, and checks')
                
                else:
                    print("Error: Unknown sheet type")
                    return
                
                ws = wb['1']
                
                jglocs = ['O16', 'O17', 'O18', 'O19', 'O20', 'O21', 'O22', 'O23', 'O24', 'O25'] # Experimental Port 1,2,3,4,5A,5C,5B,6A,6,7 --> Now Port 1,2,3,4,5,6,7,8,9,10. Quan 10/28
                try:
                    jgloc = ws[jglocs[int(re.findall(r'\d+', port)[0])-1]].value
                    #jgloc = ws[jglocs[int(port.strip('P'))]].value
                except Exception as e:
                    print(e)
                    print(f"Warning: Could not identify port # for {file}, setting jgloc = jgref")
                    jgloc = jgref

                # Above jglocs not implemented for Ryan templates (DHK)
                if jgloc is None and sheet_type == 'ryan_template':
                    print(f"Sheet type identified as {sheet_type}, referencing cell U23 for jgloc")
                    jgloc = ws['U23'].value
                elif jgloc is None:
                    print(f"Warning: jgloc could not be found, setting jgloc = jgref")
                    jgloc = jgref
                # print(theta)
                newCond = Condition(jgref, jgloc, jf, theta, port, sheet_type.split('_')[0])

                if newCond not in all_conditions:
                    all_conditions.append(newCond)
                    cond = newCond
                else:
                    cond = all_conditions[ all_conditions.index(newCond) ]

                cond.run_ID = ws['B2'].value

                # Local corrected gauge pressure can also be back-calculated from jgloc and jgatm (DHK)
               # cond.jgatm = ws['D6'].value  # Recorded experimental jgatm when data was taken, might vary slightly on different days or at different ports due to different temp or other BCs
                cond.jgatm = ws['D7'].value   # Quan 10/25 U-bend data, fixed jgatm for all the ports to avoid cofusion in later modeling 

                ws = wb['2']
                
                if sheet_type == 'ryan_template':
                    cond.area_avg_void_sheet = ws['G266'].value
                    cond.area_avg_ai_sheet = ws['J266'].value
                    
                elif sheet_type == 'adix_template' or sheet_type == 'adix_template4':
                    cond.area_avg_void_sheet = ws['G246'].value
                    cond.area_avg_ai_sheet = ws['J246'].value
                
                for phi, indices in Q1_ranges:
                    for i in indices:
                        if ws[f'{Q1_check}{i}'].value:

                            try:
                                roverR = float(ws[f'{Q1_start}{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column {Q1_start} could not be floatified. Skipping...', file=debugFID)
                                continue
                            
                            midas_output = []
                            for cell in ws[f'{Q1_start}{i}':f'{Q1_end}{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output", file=debugFID)
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)
                            
                            try:
                                cond.data[phi].update({roverR: data})
                            except KeyError:
                                cond.data.update( {phi:{}} )
                                cond.data[phi].update({1.0: deepcopy(zero_data)})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.data[phi].update({roverR: data})

                
                if pitot_sheet: # Pitot data for Q1
                    for phi, indices in Q1_ranges:
                        for i in indices:
                            if ws[f'{Q1_pitot_check}{i}'].value and ws[f'{Q1_check}{i}'].value:

                                try:
                                    roverR = float(ws[f'{Q1_pitot_start}{i}'].value)
                                except:
                                    if debug: print(f'Warning: data found in row {i} in sheet {file}, but column {Q1_start} could not be floatified. Skipping...')
                                    continue

                                pitot_output = []
                                for cell in ws[f'{Q1_pitot_start}{i}':f'{Q1_pitot_end}{i}'][0]:
                                    pitot_output.append(cell.value)

                                if len(pitot_keys2) == len( pitot_output ):
                                    pitot_data = dict( zip( pitot_keys2, pitot_output ))
                                else:
                                    if debug: print("Warning, pitot_keys2 not the same length as pitot_output", file=debugFID)
                                    pitot_data = dict( zip( pitot_keys2, pitot_output ))
                                    if debug:
                                        print("pitot_keys2 not the same length as pitot_output", file=debugFID)
                                        print(pitot_keys2, pitot_output, file=debugFID)

                                try:
                                    cond.data[phi][roverR].update(pitot_data)
                                except KeyError:
                                    cond.data[phi].update({roverR: deepcopy(zero_data)})
                                    cond.data[phi][roverR].update(pitot_data)

                for phi, indices in Q2_ranges:
                    for i in indices:
                        if ws[f'{Q2_check}{i}'].value:

                            try:
                                roverR = float(ws[f'{Q2_start}{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column {Q2_start} could not be floatified. Skipping...')
                                continue
                            
                            midas_output = []
                            for cell in ws[f'{Q2_start}{i}':f'{Q2_end}{i}'][0]:
                                midas_output.append(cell.value)
                            
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                if debug: print("Warning, tab_keys not the same length as midas_output")
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)

                            try:
                                cond.data[phi].update({roverR: data})
                            except KeyError:
                                cond.data.update( {phi:{}} )
                                cond.data[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.data[phi].update({roverR: data}) 

                if pitot_sheet: # Pitot data for Q2
                    for phi, indices in Q2_ranges:
                        for i in indices:
                            if ws[f'{Q2_pitot_check}{i}'].value:

                                try:
                                    roverR = float(ws[f'{Q2_pitot_start}{i}'].value)
                                except:
                                    if debug: print(f'Warning: data found in row {i} in sheet {file}, but column {Q2_pitot_start} could not be floatified. Skipping...')
                                    continue

                                pitot_output = []
                                for cell in ws[f'{Q2_pitot_start}{i}':f'{Q2_pitot_end}{i}'][0]:
                                    pitot_output.append(cell.value)

                                if len(pitot_keys2) == len( pitot_output ):
                                    pitot_data = dict( zip( pitot_keys2, pitot_output ))
                                else:
                                    if debug: print("Warning, pitot_keys2 not the same length as pitot_output", file=debugFID)
                                    pitot_data = dict( zip( pitot_keys2, pitot_output ))
                                    if debug:
                                        print("pitot_keys2 not the same length as pitot_output", file=debugFID)
                                        print(pitot_keys2, pitot_output, file=debugFID)

                                try:
                                    cond.data[phi][roverR].update(pitot_data)
                                except KeyError:
                                    cond.data[phi].update({roverR: deepcopy(zero_data)})
                                    cond.data[phi][roverR].update(pitot_data)

    if debug and False:
        for cond in all_conditions:
            cond.pretty_print()

    if append_to_json:
        # Read in previous data, and append this new data to that
        with open(append_to_json, 'rb') as g: # pickle works with binary data, hence rb
            data = pickle.load(g)
        
        all_conditions = data + all_conditions # In theory could produce a database with heterogeneous Condition objects, if messing with source code
    
    with open(dump_file, 'wb') as g:
        pickle.dump(all_conditions, g)
    return


def dump_data_from_tabs(dump_file = 'PITA_Database.dat', skip_dir = "") -> None:
    all_conditions = []
    
    for path, directories, files in os.walk('.'):
        

        if len(path.split('\\')) == 5:
            for file in files:
                if file.split('.')[-1] == 'tab':
                    
                    #print(path)
                    try:
                        jf = float(path.split('\\')[2].split('_')[0].strip('jf'))
                        jg = float(path.split('\\')[2].split('_')[1].strip('jg'))
                        port = path.split('\\')[3]
                        phi_val = float(path.split('\\')[4])
                        theta = float(path.split('\\')[1])

                        newCond = Condition(jg, jg, jf, theta, port)

                        if newCond not in all_conditions:
                            all_conditions.append(newCond)
                                                
                        cond = all_conditions[ all_conditions.index(newCond) ]

                        with open(os.path.join(path, file), 'r') as f:
                            midas_output = [float(i) for i in  f.readline().strip().split('\t') ]
                            if len(tab_keys) == len( midas_output ):
                                data = dict( zip( tab_keys, midas_output ))
                            else:
                                #print("Warning, tab_keys not the same length as midas_output")
                                data = dict( zip( tab_keys, midas_output ))
                                if debug:
                                    print("tab_keys not the same length as midas_output", file=debugFID)
                                    print(tab_keys, midas_output, file=debugFID)
                        roverR = float(midas_output[0])
                        if 'i' in file:
                            roverR = -roverR
                        #print(roverR)
                        try:
                            cond.data[phi_val].update({roverR: data})
                        except KeyError:
                            cond.data.update( {phi_val:{}} )
                            cond.data[phi_val].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.data[phi_val].update({roverR: data})                    
                        

                    except Exception as e:
                        if debug:
                            print(e, file=debugFID)
                            print("Continuing...", file=debugFID)

            

    if debug and True:
        print("Printing all things read from files", file=debugFID)
        for cond in all_conditions:
            cond:Condition
            if cond.theta == 90:
                cond.pretty_print()
            
                    
    
    # print("Found the following conditions:")
    # for cond in all_conditions:
    #     print(cond)

    with open(dump_file, 'wb') as g:
        pickle.dump(all_conditions, g)

    return

def loadData(data_file) -> list:
    with open(data_file, 'rb') as g: # pickle works with binary data, hence rb
        data = pickle.load(g)
    return data

def loadProbeData(data_file = 'PITA_Database.dat') -> list:
    with open(data_file, 'rb') as g: # pickle works with binary data, hence rb
        data = pickle.load(g)
    return data

def dumpData(data, dump_file = "database.dat"):
    with open(dump_file, 'wb') as g:
        pickle.dump(data, g)
    return

def extractPitotData(dump_file = 'Pitot_Database.dat', in_dir = [], require_terms = [], skip_terms = ['CFD', 'Copy']) -> None:
    Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 33)], [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)], [i for i in range(208, 233)] ]))

    path = "Z:\\TRSL\\PITA\\Data\\PitotData"

    all_conditions = []

    for file in os.listdir(path):
        if debug: print(file, file=debugFID)
        if file.split('.')[-1] == 'xlsx':
            
            # Check if the file has any skipped/required terms
            if any(term in file for term in skip_terms):
                if debug: print(f"Skipping {file}", file=debugFID)
                continue

            if any(term not in file for term in require_terms):
                if debug: print(f"Skipping {file}", file=debugFID)
                continue
            
            #if debug: print(path, file=debugFID)
            
            try:
                wb = op.load_workbook(filename=os.path.join(path, file), data_only=True)
            except:
                print(f"Error reading wb: {file}\nSkipping...")
                continue
            
            try:
                jf = float(file.split('_')[1].strip('jf'))
                jgref = float(file.split('_')[2].strip('jg'))
                port = file.split('_')[3].strip('.xlsx')
                theta = float(file.split('_')[0].strip('deg'))
            except:
                print(f'Warning: Non-standard excel file name {file}. Skipping...')
                continue

            ws = wb['1']
            jgloc = ws['U23'].value

            newCond = Condition(jgref, jgloc, jf, theta, port, 'Pitot')

            if newCond not in all_conditions:
                all_conditions.append(newCond)
                cond = newCond
            else:
                cond = all_conditions[ all_conditions.index(newCond) ]
            
            ws = wb['2']

            cond.area_avg_void_sheet = ws['G266'].value
            
            for phi, indices in Q1_ranges:
                for i in indices:
                    if ws[f'E{i}'].value:

                        try:
                            roverR = float(ws[f'A{i}'].value)
                        except:
                            if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...')
                            continue
                        
                        midas_output = []
                        for cell in ws[f'B{i}':f'BD{i}'][0]:
                            midas_output.append(cell.value)
                        
                        if len(tab_keys) == len( midas_output ):
                            data = dict( zip( pitot_keys, midas_output ))
                        else:
                            if debug: print("Warning, pitot_keys not the same length as midas_output")
                            data = dict( zip( pitot_keys, midas_output ))
                            if debug:
                                print("pitot_keys not the same length as midas_output", file=debugFID)
                                print(pitot_keys, midas_output, file=debugFID)

                        try:
                            cond.data[phi].update({roverR: data})
                        except KeyError:
                            cond.data.update( {phi:{}} )
                            cond.data[phi].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.data[phi].update({roverR: data})


    with open(dump_file, 'wb') as g:
        #print(all_conditions)
        pickle.dump(all_conditions, g)
    return

def loadPitotData(data_file = 'Pitot_Database.dat') -> list:
    with open(data_file, 'rb') as g: # pickle works with binary data, hence rb
        data = pickle.load(g)
    return data

def extractIskandraniData(dump_file='Iskandrani_Database.dat') -> None:
    try:
        wb = op.load_workbook("Z:\\TRSL\\PITA\Data\\LocalData\\Iskandrani\\Iskandrani_Local_Velocity.xlsx")
    except IOError as e:
        print(e)
        print("The path to the Iskandrani data is hardcoded, maybe the file moved?")
        return

    database = []
    for ws in wb.worksheets:
        #print(ws)
        jf = ws['B1'].value
        jg = ws['B2'].value
        new_cond = Iskandrani_Condition(jf, jg)

        new_cond.data.update({'rs': [cell[0].value for cell in ws['A4:A24']]})
        new_cond.data.update({'vf': [cell[0].value for cell in ws['B4:B24']]})
        new_cond.data.update({'alpha': [cell[0].value for cell in ws['D4:D24']]})
        new_cond.data.update({'vfp': [cell[0].value for cell in ws['C4:C24']]})
        new_cond.data.update({'bubble_freq': [cell[0].value for cell in ws['E4:E24']]})
        new_cond.data.update({'fluctuations': [cell[0].value for cell in ws['F4:F24']]})

        database.append(new_cond)

    with open(dump_file, 'wb') as g:
        pickle.dump(database, g)

    return

def loadIskandraniData(data_file = 'Iskandrani_Database.dat') -> list:
    with open(data_file, 'rb') as g: # pickle works with binary data, hence rb
        data = pickle.load(g)
    return data

def extractYangData(dump_file='Yang_Database.dat') -> list:

    try:
        wb = op.load_workbook("Z:\\TRSL\\PITA\Data\\LocalData\\Yang_et_al_mean_velocity_profiles.xlsx")
    except IOError as e:
        print(e)
        print("The path to the Yang data is hardcoded, maybe the file moved?")
        return
    
    ws = wb.active

    database = []

    for header_cell in ws['1']:

        if header_cell.value:
            jf = float(header_cell.value.split('_')[2].strip('jf='))
            jg = float(header_cell.value.split('_')[1].strip('jg='))
            phi_value = int(header_cell.value.split('_')[0])

            if database:
                if (database[-1].jf != jf) and (database[-1].jg != jg):
                    database.append( Yang_Condition(jf, jg) )
            else:
                database.append( Yang_Condition(jf, jg) )

            r_dict = {1.0: zero_data}
            for i in range(2, 12):
                r = ws.cell(row = i, column = header_cell.column).value
                vf = ws.cell(row = i, column = header_cell.column+1).value
                midas_dict = {'vf': vf}
                r_dict.update({r: midas_dict})
            database[-1].phi.update({phi_value: r_dict})      

    with open(dump_file, 'wb') as g:
        pickle.dump(database, g)

    return

def loadYangData(data_file = 'Yang_Database.dat') -> list:
    with open(data_file, 'rb') as g: # pickle works with binary data, hence rb
        data = pickle.load(g)
    return data

tab_keys = [
    'signed_roverR',
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

pitot_keys = [
    'roverR',
    'time',
    'frequency',
    'delta_p',
    'sigma_delta_p',
    'vf']

pitot_keys2 = [
    'signed_roverR',
    'roverR',
    'time',
    'frequency',
    'delta_p',
    'sigma_delta_p',
    'vf_naive',
    'vf',
    'jf',
    'vr'
]