###########################################################################################################################################
##
## A script to analyze the relationship between new parameters and the final closeness of the energy value obtained by Maestro,
## and generate new parameters based on previous data
##
###########################################################################################################################################

import sys
import os
import csv
import pandas as pd
pd.set_option('display.max_columns', None)
#import numpy as np
from pandas import DataFrame,Series
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem import TorsionFingerprints
import time
import glob
from math import sqrt
from sklearn import metrics
from sklearn.metrics import mean_absolute_error # MAE
from sklearn.metrics import mean_squared_error # MSE
from sklearn.metrics import r2_score # R square
from PIL import Image
from numpy import average, dot, linalg
from skimage.metrics import structural_similarity
from bayes_opt import BayesianOptimization
import random

start=time.time()
print("start_time:", str(start))


## A function to add hydrogen for all diamond energy direct conformers output
def diamond_addhs(file_order, energy_file_order):
    ## file_order here means molecule file order, energy_file_order means the order of current using energy txt file
    ## The direct conformers output path
    ## Need to read the energy_list stored the lowest 25 energy conformers number
    if os.path.exists(str(file_order)+"_"+str(energy_file_order)+"_energy_list.txt"):
        with open(str(file_order)+"_"+str(energy_file_order)+"_energy_list.txt", "r+") as f:
            for line in f.readlines():
                line=line.strip()
                print("Current_conf_num:", line)
                if os.path.exists(str(file_order)+"_"+str(line)+"_"+str(energy_file_order)+"_.sdf"):
                    suppl_mol=Chem.SDMolSupplier(str(file_order)+"_"+str(line)+"_"+str(energy_file_order)+"_.sdf")
                    mol=suppl_mol[0]
                    mol_addhs=Chem.AddHs(mol, addCoords=True)
                    mol_addhs_file=Chem.MolToMolFile(mol_addhs, "H"+str(file_order)+"_"+str(line)+"_"+str(energy_file_order)+"_.sdf")



##for each .mol2  file  in dataset, convert them into .mae format that could be used in Maestro
def sdf_to_mae(file_order, energy_file_order):
    # Pattern to match files of interest, assuming `file_order` and `energy_file_order` are known and fixed
    pattern = f"H{file_order}_*_{energy_file_order}_.sdf"

    # Use glob to find all .sdf files matching the pattern
    sdf_files = glob.glob(pattern)

    for sdf_file in sdf_files:
        # Extract the conf_num part from the file name
        parts = sdf_file.split('_')
        conf_num = parts[1]

        # Construct the script file name
        script_file_name = f"sdf_to_mae_{file_order}_{conf_num}_{energy_file_order}.sh"

        with open(script_file_name, "w+") as f:
            f.write("#!/bin/bash\n")
            f.write(
                f"/shared/shared/schrodinger/2022-1/utilities/structconvert -isd {sdf_file} -omae mini_{file_order}_{conf_num}_{energy_file_order}.mae\n")

        # Execute the shell script
        command = f"bash {script_file_name}"
        os.system(command)
        print(f"Processed {sdf_file}")


## .com file for minimization task
def mini_file_write(file_order, energy_file_order):
    # Pattern to match files of interest, assuming file_order and energy_file_order are known and fixed
    pattern = f"mini_{file_order}_*_{energy_file_order}.mae"
    # Use glob to find all matching .mae files
    for mae_file in glob.glob(pattern):
        # Extract the conf_num from the file name
        parts = mae_file.split('_')
        conf_num = parts[1]  # Assuming the file name format is "mini_fileOrder_confNum_energyFileOrder.mae"
        com_file_name = f"mini_{file_order}_{conf_num}_{energy_file_order}.com"
        with open(com_file_name, "w") as f4:
            f4.write(f"{mae_file}\n")
            f4.write(f"mini_{file_order}_{conf_num}_{energy_file_order}-out.maegz\n")
            f4.writelines([
                " MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n",
                " DEBG      55      0      0      0     0.0000     0.0000     0.0000     0.0000\n",
                " DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n",
                " FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n",
                " SOLV       3      1      0      0     0.0000     0.0000     0.0000     0.0000\n",
                " EXNB       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n",
                " BDCO       0      0      0      0    89.4427 99999.0000     0.0000     0.0000\n",
                " CRMS       0      0      0      0     0.0000     0.5000     0.0000     0.0000\n",
                " BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n",
                " READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n",
                " CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000\n",
                " MINI       1      0   2500      0     0.0000     0.0000     0.0000     0.0000\n",
                " END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n",
            ])



##for each mol wile  a structure file .mae  and a command file .com
##write the bash file for them
def write_CSbash_file(file_order, energy_file_order):
    pattern = f"mini_{file_order}_*_{energy_file_order}.com"
    # Use glob to find all matching .mae files
    for com_file in glob.glob(pattern):
        # Extract the conf_num from the file name
        parts = com_file.split('_')
        conf_num = parts[1]  # Assuming the file name format is "mini_fileOrder_confNum_energyFileOrder.com"
        with open("mini_"+str(file_order)+"_"+str(conf_num)+"_"+str(energy_file_order)+".sh","w+") as f2:
            f2.write("#!/bin/bash")
            f2.write("\n")
            f2.write("/shared/shared/schrodinger/2022-1/bmin -WAIT mini_"+str(file_order)+"_"+str(conf_num)+"_"+str(energy_file_order))
            f2.write("\n")
            #f2.write("/shared/shared/schrodinger/2019-1/utilities/structconvert -imae mini_"+str(mol_order)+"-out.maegz -osd mini_"+str(mol_order)+"_result.sdf")
        f2.close()


##for each bash file for minimization in MacroModel, run them one by one
def run_CS_onebyone(file_order, energy_file_order):
    pattern = f"mini_{file_order}_*_{energy_file_order}.sh"
    # Use glob to find all matching .mae files
    for sh_file in glob.glob(pattern):
        # Extract the conf_num from the file name
        parts = sh_file.split('_')
        conf_num = parts[1]  # Assuming the file name format is "mini_fileOrder_confNum_energyFileOrder.sh"
        ## The current working path
        if os.path.exists("mini_"+str(file_order)+"_"+str(conf_num)+"_"+str(energy_file_order)+".sh"):
            command="bash mini_"+str(file_order)+"_"+str(conf_num)+"_"+str(energy_file_order)+".sh"
            os.system(command)


##for each output .maegz file, convert them into .sdf format for further comparison with the results of Diamond_energy program
def convert2sdf_file(file_order, energy_file_order):
    pattern = f"mini_{file_order}_*_{energy_file_order}-out.maegz"
    # Use glob to find all matching .mae files
    for maegz_file in glob.glob(pattern):
        # Extract the conf_num from the file name
        parts = maegz_file.split('_')
        conf_num = parts[1]  # Assuming the file name format is "mini_fileOrder_confNum_energyFileOrder.sh"
        ## the current working path
        with open("convert2sdf_"+str(file_order)+"_"+str(conf_num)+"_"+str(energy_file_order)+".sh","w+") as f3:
            f3.write("#!/bin/bash")
            f3.write("\n")
            f3.write("/shared/shared/schrodinger/2022-1/utilities/structconvert -imae mini_"+str(file_order)+"_"+str(conf_num)+"_"+str(energy_file_order)+"-out.maegz -osd mini_result_"+str(file_order)+"_"+str(conf_num)+"_"+str(energy_file_order)+".sdf")
        f3.close()
        command="bash convert2sdf_"+str(file_order)+"_"+str(conf_num)+"_"+str(energy_file_order)+".sh"
        os.system(command)


## read the data from previous .csv results
## choose where the data needed from their index
#df=pd.read_csv("test.csv", sep=',', usecols=['para1', 'para2', 'para3', 'para4', 'para5', 'para6', 'para7', 'para8', 'para9', 'para10', 'para11', 'para12','para13', 'para14', 'para15', 'R2'])
#print(df)
#print('Shape:',df.shape)
## get the initial parameters
current_row=0
best_value_row=0
current_para_order=0
total_para_num=15
standard_row=best_value_row


# Define optimizer to optimise parameters
# Set search space/bounds for parameters
pbounds={
    'para1': (-3, -1),
    'para2': (-3, -1),
    'para3': (2, 6),
    'para4': (-8, -4),
    'para5': (-4, -1),
    'para6': (8, 12),
    'para7': (-2, -1),
    'para8': (-3, -2),
    'para9': (1, 3),
    'para10': (-2, -1),
    'para11': (-4, -2),
    'para12': (-15, -10),
    'para13': (-9, -7),
    'para14': (-5, -3),
    'para15': (-20, -15),
    'para16': (-27, -24),
    'para17': (-32, -28),
    'para18': (-22, -18),
    'para19': (-23, -21)}


# This is the function wasnt to be minimised
global process
process=0
def objective(para1, para2, para3, para4, para5, para6, para7, para8, para9, para10, para11, para12, para13, para14, para15):#, para16, para17, para18, para19):
    global process
    file_order=process
    ## in the setting path, generate new parameters set and write in an energy file
    with open(str(file_order)+"data.txt", "a+") as f:
        n=1
        while n < 16:
            if n==1:
                energy_value=round(para1)
            elif n==2:
                energy_value=round(para2)
            elif n==3:
                energy_value=round(para3)
            elif n==4:
                energy_value=round(para4)
            elif n==5:
                energy_value=round(para5)
            elif n==6:
                energy_value=round(para6)
            elif n==7:
                energy_value=round(para7)
            elif n==8:
                energy_value=round(para8)
            elif n==9:
                energy_value=round(para9)
            elif n==10:
                energy_value=round(para10)
            elif n==11:
                energy_value=round(para11)
            elif n==12:
                energy_value=round(para12)
            elif n==13:
                energy_value=round(para13)
            elif n==14:
                energy_value=round(para14)
            elif n==15:
                energy_value=round(para15)
            elif n==16:
                energy_value=round(para16)
            elif n==17:
                energy_value=round(para17)
            elif n==18:
                energy_value=round(para18)
            elif n==19:
                energy_value=round(para19)
            #elif n==xx:
            #    energy_value=round(paraxx)
            print("current_energy_value:", energy_value)
            f.write(str(energy_value))
            f.write("\n")
            n+=1
    f.close()
    ## now to use this new parameter set to do conformaitonal searching in diamond_energy script
    mol_order=1
    diamond_energy_collect=[]
    maestro_energy_collect=[]
    diamond_energy_lists=[]
    maestro_energy_lists=[]
    for smile in smile_list:
        mol_input = smile.strip()
        #real_mol = Chem.MolToInchi(Chem.MolFromSmiles(smile_input))
        ## test them by using current energy parameter set file,
        ## the ouput files generated by diamond energy will then be names with,
        ## parameter set file order_conformers generated order_.sdf
        command = "python diamond_energy.py " + "\"" + mol_input + "\"" + " " + str(mol_order) + " " + str(file_order)
        os.system(command)
        ## after getting all the direct results from diamond_energy script
        ## add hydrogens in the diamond_energy output .sdf molecules
        diamond_addhs(mol_order, file_order)
        ## convert adding hydrogens files into .mea format for minimization in maestro
        sdf_to_mae(mol_order, file_order)
        ## write .com task file
        mini_file_write(mol_order, file_order)
        ## write bash file for commiting task
        write_CSbash_file(mol_order, file_order)
        ## do the minimization task
        run_CS_onebyone(mol_order, file_order)
        ## convert the minimization result into .sdf file for further align comparison
        convert2sdf_file(mol_order, file_order)
        #############################################################################################################
        ##
        ## now going to collect the data got from diamond_energy for new round's analysis
        ##
        ##
        #############################################################################################################
        print("\n")
        print("Below is diamond's result:")
        print("\n")

        ## To speed up the comparison, now cut the total energy numbers to 25
        def process_energy_files(mol_order, file_order):
            # Create a pattern that matches all relevant log files
            pattern = f"mini_{mol_order}_*_{file_order}.log"
            # Use glob to find all files matching the pattern
            log_files = glob.glob(pattern)
            diamond_energy = []
            # Process each file found by glob
            for log_file in log_files:
                with open(log_file, "r") as f:
                    num_TE = 0
                    for line in f:
                        if "T.E. for cross-checking" in line:
                            num_TE += 1
                            if num_TE == 2:
                                current_energy = float(line.split(":")[1].split("kJ")[0])
                                diamond_energy.append(current_energy)
                                break  # Exit the loop after the second occurrence

            return diamond_energy
        diamond_energy = process_energy_files(mol_order, file_order)
        ## if the results of diamond is none, means that the parameter makes it not working.
        ## in this situation to avoid the comparison stop for none data from diamond\
        ## to create any artificial list such as [0.1, 0.1, 0.1, 0.1, 0.1] as inout
        if not diamond_energy:
            diamond_energy=[0.1, 0.1, 0.1, 0.1, 0.1]
        else:
            ## as the energy might not sorted in order so need to sort the energy list got
            diamond_energy.sort()
            diamond_energy=diamond_energy[:5]
            lowest_dia_energy=min(diamond_energy)
            print("lowest_dia_energy:", lowest_dia_energy)
            number_of_low_energy_conformer_dia=len(diamond_energy)
            print("number_of_low_energy_conformer_dia:", number_of_low_energy_conformer_dia)
            ## calculate the second lowest energy number
            second_num_dia = 0
            for value in diamond_energy:
                # print("value:", value)
                if value <= (lowest_dia_energy + 1):
                    second_num_dia += 1
            print("second_num_dia:", second_num_dia)

        ##############################################################################################################
        ##
        ## now going to collect the data got from Maestro as benchmark
        ##
        ##############################################################################################################

        print("\n")
        print("Below is maestro's result:")
        print("\n")

        with open("/home/mw742/Meastro/mmod_csearch_"+str(mol_order)+"/mmod_csearch_"+str(mol_order)+".log", 'r') as f:
            f_txt = f.read().splitlines()
            # print(type(f_txt))
            # print("f_txt:", f_txt)
            final_report_index = f_txt.index("Final report:")
            # print("final_report_index:", final_report_index)
            num_of_uni_stru = int(f_txt[final_report_index + 1].split("unique")[0])
            print("num_of_uni_stru:", num_of_uni_stru)

            global_mini_line = [i for i in f_txt if "Conformation       1" in i]
            global_mini_index = f_txt.index(global_mini_line[0])
            # print("global_mini_index:", global_mini_index)
            num_of_global_mini = int(f_txt[global_mini_index].split("found")[1].split("times")[0])
            global_mini_energy = f_txt[global_mini_index].split("(")[1].split(")")[0]
            print("freq_of_finding_global_mini:", num_of_global_mini)
            print("global_mini_energy:", global_mini_energy)

            ## count the structures which are within 1 kJ/mol of the global minimum, what is the smallest number of times any of them has been found
            num_within_mini_range_index = final_report_index + 3
            num_within_mini_range = int(f_txt[num_within_mini_range_index].split("confs")[0].split("Found")[1])
            a = 0
            current_mini_num = 0
            current_index = global_mini_index
            smallest_mini_num = num_of_global_mini
            while a < num_within_mini_range:
                current_mini_num = int(f_txt[current_index].split("found")[1].split("times")[0])
                if current_mini_num < smallest_mini_num:
                    smallest_mini_num = current_mini_num
                else:
                    pass
                current_index += 1
                a += 1
            print("smallest_mini_num:", smallest_mini_num)

            ## calculate the time to run a search
            time_line = [i for i in f_txt if "Time in Monte Carlo generation loop:" in i]
            time_index = f_txt.index(time_line[0])
            conformer_generate_time = float(f_txt[time_index].split(":")[1].split("CPU")[0])
            energy_minimization_time = float(f_txt[time_index + 1].split(":")[1].split("CPU")[0])
            total_time = conformer_generate_time + energy_minimization_time
            print("total_time:", total_time)

            ##collect all the energy value in a list for further comparison
            benchmark_energy = []
            for line in f_txt:
                if line.find(") was found") > -1:
                    current_energy_value = float(line.split("kJ/mol")[0].split("(")[1])
                    # print("current_energy_value:", current_energy_value)
                    benchmark_energy.append(current_energy_value)
            print("number_of_energy_value:", len(benchmark_energy))
            ## to speed up the comparison, now cut the total energy numbers to 20
            benchmark_energy.sort()
            lowest_benchmark_energy = min(benchmark_energy)
            number_of_benchmark_conformer = len(benchmark_energy)
            number_of_benchmark_second = smallest_mini_num
            f.close()
        ## here to compare the results of two methods
        ## make the length of the comparison limit to below 25
        len_maestro = len(benchmark_energy)
        len_diamond = len(diamond_energy)
        if len_diamond >= len_maestro:
            num_compare = len_maestro
            diamond_energy_final = diamond_energy[:num_compare]
            benchmark_energy_final = benchmark_energy
        else:
            num_compare = len_diamond
            benchmark_energy_final = benchmark_energy[:num_compare]
            diamond_energy_final = diamond_energy
        diamond_energy_collect.extend(diamond_energy_final)
        maestro_energy_collect.extend(benchmark_energy_final)
        diamond_energy_lists.extend([diamond_energy_final])
        maestro_energy_lists.extend([benchmark_energy_final])
        mol_order+=1

    # Convert to np.array since this is way faster than processing a list of lists
    # Formula for mse
    # Can also use sklearn package if you want
    # Can change metric to minimise by just changing the line below and the respective variables
    # mse = ((a-b)**2).mean(axis=1)
    # maestro=np.array(maestro_energy_collect)
    # diamond=np.array(diamond_energy_collect)
    R2=metrics.r2_score(maestro_energy_collect, diamond_energy_collect)
    process+=1
    return R2




##############################################################################################################
##
## Starting point parameter set value and range
##
##############################################################################################################



## read the test molecules dataset from prepared Test_list.txt
global smile_list
smile_list=[]
with open("Test_list.txt", "r") as f:
    for line in f.readlines():
        smile_list.append(line)
f.close()




#############################################################################################################################
##
## The automatically analysis and re-generate parameters set to loop through para-pcc surface starts from here
##
#############################################################################################################################



# Run optimiser to find best parameter set to use
optimizer=BayesianOptimization(
    f=objective,
    pbounds=pbounds,
    verbose=2,
    random_state=2)


best_params=optimizer.maximize(init_points=15, n_iter=100)
print("best_params:", optimizer.max)


# This should give you a list of the best parameter set to use print(best_params)
# To test best_params
# R2=objective(optimizer.max)


end=time.time()
print("end_time:", str(end))
print("running_time", str(end-start))
