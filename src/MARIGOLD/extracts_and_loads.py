
from .config import *
from .Condition import Condition
from .Iskandrani_Condition import Iskandrani_Condition
from .Yang_Condition import Yang_Condition

import re

def extractProbeData(dump_file = 'database.dat', in_dir = [], require_terms = None, skip_terms = ['CFD', 'Copy']) -> None:
    all_conditions = []

    ### PITA Data ###

    path = 'Z:\\TRSL\\PITA\\Data\\LocalData\\spreadsheets\\PITA'
    #for path, directories, files in os.walk('./spreadsheets'):
    #for file in os.listdir(path):
    if True:

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
                    jgP3 = float(file.split('_')[2].strip('jg'))
                    port = file.split('_')[3].strip('.xlsx')
                    theta = float(file.split('_')[0].strip('deg'))
                except:
                    print(f'Warning: Non-standard excel file name {file}. Skipping...')
                    continue

                ws = wb['1']
                jgloc = ws['U23'].value

                newCond = Condition(jgP3, jgloc, jf, theta, port, 'Ryan')

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
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...', file=debugFID)
                                continue
                            
                            midas_output = []
                            for cell in ws[f'B{i}':f'BD{i}'][0]:
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
                                cond.phi[phi].update({roverR: data})
                            except KeyError:
                                cond.phi.update( {phi:{}} )
                                cond.phi[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.phi[phi].update({roverR: data})


                for phi, indices in Q2_ranges:
                    for i in indices:
                        if ws[f'DC{i}'].value:

                            try:
                                roverR = float(ws[f'CO{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...')
                                continue
                            
                            midas_output = []
                            for cell in ws[f'CP{i}':f'ER{i}'][0]:
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
                                cond.phi[phi].update({roverR: data})
                            except KeyError:
                                cond.phi.update( {phi:{}} )
                                cond.phi[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.phi[phi].update({roverR: data})

    
    ### 38.1 mm data (Talley) ###
    path = 'Z:\\TRSL\\PITA\\Data\\LocalData\\spreadsheets\\38.1mm'

    if True:

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
                    jgP3 = float(file.split('_')[2].strip('jg'))
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
                

                newCond = Condition(jgP3, jgloc, jf, theta, port, 'Talley')

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
                            for cell in ws[f'B{i}':f'BD{i}'][0]:
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
                            cond.phi[phi].update({roverR: data})
                        except KeyError:
                            cond.phi.update( {phi:{}} )
                            cond.phi[phi].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.phi[phi].update({roverR: data})

                #print('total rows read:', i)

    ### 101.6mm Data ###

    path = 'Z:\\TRSL\\PITA\\Data\\LocalData\\spreadsheets\\101.6mm'
    #for path, directories, files in os.walk('./spreadsheets'):
    #for file in os.listdir(path):
    if True:

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
                    jgP3 = float(file.split('_')[2].strip('jg'))
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
                

                newCond = Condition(jgP3, jgloc, jf, theta, port, 'Kong')

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
                            for cell in ws[f'B{i}':f'BD{i}'][0]:
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
                                cond.phi[phi].update({roverR: data})
                            except KeyError:
                                cond.phi.update( {phi:{}} )
                                cond.phi[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.phi[phi].update({roverR: data})


                for phi, indices in Q2_ranges:
                    for i in indices:
                        if ws[f'DC{i}'].value:

                            try:
                                roverR = float(ws[f'CN{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column A could not be floatified. Skipping...')
                                continue
                            
                            midas_output = []
                            for cell in ws[f'CO{i}':f'EQ{i}'][0]:
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
                                cond.phi[phi].update({roverR: data})
                            except KeyError:
                                cond.phi.update( {phi:{}} )
                                cond.phi[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.phi[phi].update({roverR: data})
                
    
    ### Yadav Data ###

    path = 'Z:\\TRSL\\PITA\\Data\\LocalData\\spreadsheets\\50.8mm'
    if True:
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
                    jgP3 = float(file.split('_')[2].strip('jg'))
                    port = file.split('_')[3].strip('.xlsx')
                    theta = float(file.split('_')[0].strip('deg'))
                except:
                    print(f'Warning: Non-standard excel file name {file}. Skipping...')
                    continue

                ws = wb['1']

                # jglocs = ['O16', 'O17', 'O18', 'O19', 'O20']
                # jgloc = ws[jglocs[int(port.strip('P'))-1]].value
                jgloc = jgP3
                
                newCond = Condition(jgP3, jgloc, jf, theta, port, 'Yadav')

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
                            for cell in ws[f'B{i}':f'BD{i}'][0]:
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
                                cond.phi[phi].update({roverR: data})
                            except KeyError:
                                cond.phi.update( {phi:{}} )
                                cond.phi[phi].update({1.0: zero_data})
                                #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                                cond.phi[phi].update({roverR: data})


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
                            cond.phi[phi].update({roverR: data})
                        except KeyError:
                            cond.phi.update( {phi:{}} )
                            cond.phi[phi].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.phi[phi].update({roverR: data}) 

    
    if debug and False:
        for cond in all_conditions:
            cond.pretty_print()
    
    with open(dump_file, 'wb') as g:
        pickle.dump(all_conditions, g)
    return

def extractLocalDataFromdDir(path, dump_file = 'database.dat', in_dir = [], require_terms = ['jf'], 
                            skip_terms = ['CFD', 'Copy'], sheet_type = 'adix_template', append_to_json = None,
                            custom_ranges = None, custom_starts_and_ends = None
                            ) -> None:
    
    """Function for getting all local data from spreadhseets in a directory
    
    
    """
    all_conditions = []

    for file in os.listdir(path):
    
        if debug: print(file, file=debugFID)
        if file.split('.')[-1] == 'xlsx' or file.split('.')[-1] == 'xlsm':

            if sheet_type == 'infer':
                if file.split('.')[-1] == 'xlsm' or 'adix' in file:
                    sheet_type = 'adix_template4'
                elif '60deg' in file or '30deg' in file or '80deg' in file or 'Ryan' in file:
                    sheet_type = 'ryan_template'
                elif 'P4' in file or 'P5' in file or 'P6' in file or 'P7' in file:
                    sheet_type = 'adix_template' # Might rename this Quan template
                else:
                    sheet_type = 'adix_template'

            if sheet_type == 'ryan_template':
                Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 33)], [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)], [i for i in range(208, 233)] ]))
                Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(57, 82)], [i for i in range(108, 133)], [i for i in range(157, 182)] ]))
                Q2_start = 'CQ'
                Q2_end = 'ES'
                Q1_start = 'B'
                Q1_end = 'BD'
            
            elif sheet_type == 'adix_template':
                Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 29)], [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)], [i for i in range(192, 213)] ]))
                Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)] ]))
                Q2_start = 'CQ'
                Q2_end = 'ES'
                Q1_start = 'B'
                Q1_end = 'BD'
            
            elif sheet_type == 'adix_template4':
                Q1_ranges = list(zip([90, 67.5, 45, 22.5, 0], [ [i for i in range(8, 29)], [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)], [i for i in range(192, 213)] ]))
                Q2_ranges = list(zip([112.5, 135, 157.5], [ [i for i in range(53, 74)], [i for i in range(100, 121)], [i for i in range(145, 166)] ]))
                Q2_start = 'CR'
                Q2_end = 'ET'
                Q1_start = 'B'
                Q1_end = 'BD'
                pitot_start = 'CG'
                pitot_end = 'CO'

            elif sheet_type == 'custom' or sheet_type == 'Custom':
                Q1_ranges = custom_ranges[0]
                Q2_ranges = custom_ranges[1]

                Q1_start = custom_starts_and_ends[0]
                Q1_end =   custom_starts_and_ends[1]
                Q2_start = custom_starts_and_ends[2]
                Q2_end =   custom_starts_and_ends[3]

                pitot_start = custom_starts_and_ends[4]
                pitot_end = custom_starts_and_ends[5]
            
            else:
                print("Error: Unknown sheet type")
                return
            
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
                port = file.split('_')[3].strip('.xlsx').strip('.xlsm')
                theta = float(file.split('_')[0].strip('deg'))
            except:
                print(f'Warning: Non-standard excel file name {file}. Skipping...')
                continue

            ws = wb['1']

            jglocs = ['O16', 'O17', 'O18', 'O19', 'O20']
            try:
                jgloc = ws[jglocs[int(re.findall(r'\d+', port)[0])]].value
                #jgloc = ws[jglocs[int(port.strip('P'))]].value
            except Exception as e:
                print(e)
                print(f"Warning: Could not identify port # for {file}, setting jgloc = jgref")
                jgloc = jgref

            newCond = Condition(jgref, jgloc, jf, theta, port, sheet_type.split('_')[0])

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
                            cond.phi[phi].update({roverR: data})
                        except KeyError:
                            cond.phi.update( {phi:{}} )
                            cond.phi[phi].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.phi[phi].update({roverR: data})

            
            if sheet_type == 'adix_template4':
                for phi, indices in Q1_ranges:
                    for i in indices:
                        if ws[f'CJ{i}'].value:

                            try:
                                roverR = float(ws[f'CF{i}'].value)
                            except:
                                if debug: print(f'Warning: data found in row {i} in sheet {file}, but column {Q1_start} could not be floatified. Skipping...')
                                continue

                            pitot_output = []
                            for cell in ws[f'{pitot_start}{i}':f'{pitot_end}{i}'][0]:
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
                                cond.phi[phi][roverR].update(pitot_data)
                            except KeyError:
                                cond.phi[phi].update({roverR: zero_data})
                                cond.phi[phi].update({roverR: pitot_data})

            for phi, indices in Q2_ranges:
                for i in indices:
                    if ws[f'DC{i}'].value:

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
                            cond.phi[phi].update({roverR: data})
                        except KeyError:
                            cond.phi.update( {phi:{}} )
                            cond.phi[phi].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.phi[phi].update({roverR: data}) 
    
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
                            cond.phi[phi_val].update({roverR: data})
                        except KeyError:
                            cond.phi.update( {phi_val:{}} )
                            cond.phi[phi_val].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.phi[phi_val].update({roverR: data})                    
                        

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
                jgP3 = float(file.split('_')[2].strip('jg'))
                port = file.split('_')[3].strip('.xlsx')
                theta = float(file.split('_')[0].strip('deg'))
            except:
                print(f'Warning: Non-standard excel file name {file}. Skipping...')
                continue

            ws = wb['1']
            jgloc = ws['U23'].value

            newCond = Condition(jgP3, jgloc, jf, theta, port, 'Pitot')

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
                            cond.phi[phi].update({roverR: data})
                        except KeyError:
                            cond.phi.update( {phi:{}} )
                            cond.phi[phi].update({1.0: zero_data})
                            #cond.phi[phi_val].update({0.0: zero_data}) # Cuz I'm paranoid
                            cond.phi[phi].update({roverR: data})


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
        wb = op.load_workbook("Z:\\TRSL\\PITA\Data\\LocalData\\Iskandrani_Local_Velocity.xlsx")
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