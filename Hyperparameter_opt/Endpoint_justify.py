###########################################################################################################################################
##
## A script to compare the closeness of the results between DiamondEenergyII script and MacroModel LowMode search method
##
###########################################################################################################################################


import os
import math
import csv
import pandas as pd
import sys
import rdkit
pd.set_option('display.max_columns', None)
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem import TorsionFingerprints

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
from rdkit.Chem import rdMolTransforms
import time
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

## A function to turn one conformer into mol format
def ConfToMol(mol, input):
    conf=input
    new_mol=Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(Chem.Conformer(conf))
    return new_mol


## A function to add hydrogen for all diamond energy direct conformers output
def diamond_addhs(file_order, energy_file_order):
    ## file_order here means molecule file order, energy_file_order means the order of current using energy txt file
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

##############################################################################################################
##
## Start to read the results and do the comparison
##
##############################################################################################################
column=["Global_min_diamond", "Global_min_maestro", "Boltzman_ratio_diamond", "Boltzman_ratio_maestro"]
df=pd.DataFrame(columns=column)
para_dict={"Global_min_diamond":[0],
            "Global_min_maestro":[0],
           "Boltzman_ratio_diamond":[0],
           "Boltzman_ratio_maestro":[0]}


## using the best parameter set to test the diamond_energy script's performance and collect all its
## conformational searching results
## read the test molecules dataset from prepared Test_list.txt
global smile_list
smile_list=[]
with open("Test_list.txt", "r") as f:
    for line in f.readlines():
        smile_list.append(line)
f.close()

## start to read the results data

file_order=0
mol_order=1
for smile in smile_list:
    smile_input=smile.strip()
    #real_mol=Chem.MolToInchi(Chem.MolFromSmiles(smile_input))
    real_mol=smile_input
    ## test them by using current energy parameter set file,
    ## the ouput files generated by diamond energy will then be names with,
    ## parameter set file order_conformers generated order_.sdf
    command="python diamond_energy.py " + "\"" + real_mol + "\"" + " " + str(mol_order) + " " + str(file_order)
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
    #############################################################################################################
    print("\n")
    print("Below is diamond's result:")
    print("\n")

    diamond_energy = []
    print("ATTENTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    with open(str(mol_order) + "_" + str(file_order) + "_energy_list.txt", "r+") as f:
        for number in f.readlines():
            number = number.strip()
            print("Current_conf_num:", number)
            m=0
            while m < 13:
                if os.path.exists("mini_" + str(mol_order) + "_" + str(number) + "_" + str(file_order)+"_"+str(m)+".log"):
                    # print("test")
                    with open("mini_" + str(mol_order) + "_" + str(number) + "_" + str(file_order)+"_"+str(m)+".log","r") as f:
                        num_TE = 0
                        for line in f.readlines():
                            if line.find("T.E. for cross-checking") > -1:
                                num_TE += 1
                                if num_TE == 2:
                                    current_energy=float(line.split(":")[1].split("kJ")[0])
                                    print(current_energy)
                                    diamond_energy.append(current_energy)
                                    break
                else:
                    pass
                m+=1
    f.close()

    ## if the results of diamond is none, means that the parameter makes it not work.
    ## in this situation to avoid the comparison stop for none data from diamond\
    ## to create any artificial list such as [0.1, 0.1, 0.1, 0.1, 0.1] as input
    if not diamond_energy:
        diamond_energy=[0.1, 0.1, 0.1, 0.1, 0.1]
    else:
        ## as the energy might not sorted in order so need to sort the energy list got
        diamond_energy.sort()
    print("diamond_read:", diamond_energy)

    ##############################################################################################################
    ##
    ## now going to collect the data got from Maestro as benchmark
    ##
    ##############################################################################################################

    print("\n")
    print("Below is maestro's result:")
    print("\n")

    with open("/home/mw742/Maestro/mmod_csearch_"+str(mol_order)+"/mmod_csearch_"+str(mol_order)+".log", 'r') as f:
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
        benchmark_energy=[]
        for line in f_txt:
            if line.find(") was found") > -1:
                current_energy_value=float(line.split("kJ/mol")[0].split("(")[1])
                print("current_energy_value:", current_energy_value)
                benchmark_energy.append(current_energy_value)
        benchmark_energy.sort()
        print("maestro_data:", benchmark_energy)
        f.close()


    ## to read the energy of diamond and maestro from results
    result_diamond_original=diamond_energy
    result_maestro_original=benchmark_energy
    ## to convert the energy unit from 1 kJ/mol to its equal 0.238846 kcal/mol
    scale_factor=0.238846
    result_diamond=[]
    for energy in result_diamond_original:
        current_energy=float(energy)*scale_factor
        result_diamond.append(current_energy)
    result_maestro=[]
    for energy in result_maestro_original:
        current_energy=float(energy)*scale_factor
        result_maestro.append(current_energy)
    ## only count the unique energy values
    result_diamond_unique=list(set(result_diamond))
    result_maestro_unique=list(set(result_maestro))
    diamond_min=min(result_diamond_unique)
    maestro_min=min(result_maestro_unique)


    ## here start to calculate the compare score for two methods
    ## to collect the data from both diamond and maestro to approach possible ground truth first
    possible_ground_truth=list(set(result_maestro_unique+result_diamond_unique))
    Energy_ref=min(possible_ground_truth)
    Q_coefficient=[]
    Q_coefficient_diamond=[]
    Q_coefficient_maestro=[]
    for each_energy in possible_ground_truth:
        Q_score=math.exp(-(each_energy-Energy_ref)/(298.15*0.00198872*0.238846))
        Q_coefficient.append(Q_score)
        if each_energy in result_diamond_unique:
            Q_coefficient_diamond.append(Q_score)
        if each_energy in result_maestro_unique:
            Q_coefficient_maestro.append(Q_score)
    Q_total=sum(Q_coefficient)
    Q_diamond_total=sum(Q_coefficient_diamond)
    Q_maestro_total=sum(Q_coefficient_maestro)
    diamond_score=Q_diamond_total/Q_total
    maestro_score=Q_maestro_total/Q_total
    df = df.append({"Global_min_diamond": diamond_min,
         "Global_min_maestro": maestro_min,
         "Boltzman_ratio_diamond": diamond_score,
         "Boltzman_ratio_maestro": maestro_score}, ignore_index=True)

    mol_order+=1


print("final_df", df)
df.to_csv("results.csv")