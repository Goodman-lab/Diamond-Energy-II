import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import Draw
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem import TorsionFingerprints
from rdkit.Chem.Draw import rdDepictor
from rdkit import rdBase
import os
import sys
import glob
import re

def count_files_in_directory(directory, pattern, exclude_substring=None):
    files = glob.glob(os.path.join(directory, pattern))
    if exclude_substring:
        files = [f for f in files if exclude_substring not in f]
    return len(files)


def list_files_in_directory(directory, extension, exclude_substring=None):
    if exclude_substring:
        return [f for f in os.listdir(directory) if f.endswith(extension) and exclude_substring not in f]
    else:
        return [f for f in os.listdir(directory) if f.endswith(extension)]


def list_files_in_directory(directory, pattern):
    return [f for f in os.listdir(directory) if re.match(pattern, f)]




# The molecular name in extension work test set
# One can change the following list to target test set
name=["e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9",
      "e10", "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19",
      "e20", "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29",
      "e30", "e31", "e32", "e33", "e34", "e35", "e36", "e37"]
total_number=len(name)


# The comparing experiment names of Macromodel with algorithm and force filed
Maestro_test=["MCMM_OPLS4", "MCMM_MM3", "SPMC_OPLS4", "SPMC_MM3",
              "LMCS_OPLS4", "LMCS_MM3", "MCLMCS_OPLS4", "MCLMCS_MM3"]
total_method=len(Maestro_test)





# Generate 3D structure from InChi
# For molecules tested need to run maestro calculation,
# Here are the lines to turn them into .sdf format for MacroModel running input
def inchi_to_sdf():
    with open("Test_list.txt", "r") as f:
        #ßprint(f.read().split())
        inchi_list=[]
        for line in f.read().split():
            print(line)
            inchi_list.append(line)
        mol_totalnum=len(inchi_list)
        print("mollist_length", mol_totalnum)
        inchi_n=0
        for mol_inchi in inchi_list:
            mol=Chem.MolFromInchi(mol_inchi)
            rdkit.Chem.Draw.MolToFile(mol, name[inchi_n] + "_graph.png", size=(300, 300))
            mol_maestro=Chem.AddHs(mol, addCoords=True)
            #AllChem.EmbedMolecule(mol_maestro)
            mol_maetsro_file=Chem.MolToMolFile(mol_maestro, name[inchi_n]+"_maestro.sdf")
            inchi_n+=1


# Generate input mol file for running MacroModel
# Convert SDF to MAE
def run_command(command):
    os.system(command)

#run_command("bash sdftomae.sh")



# Generate command file for running MacroModel
# to run MacroModel from command line for batch running, instead of using its graphic interaction software Maestro by hand
# We need a molecular structure file in .mae format, which can be converted using structconvert.py in Schrodinger/version/utilities/ directory
# and a command file in .com format, which includes all the tasks (operation codes and instructions can be found in macromodel_reference_manual)
# function for generating MacroModel command .com file for corresponding molecular structure .mae file

def generate_comfile(method_name):
    if method_name=="MCMM_OPLS4":
        # MCMM_OPLS4 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_MCMM_OPLS4.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_MCMM_OPLS4-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" MCMM    2000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
                f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
                f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_MCMM_OPLS4".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="MCMM_MM3":
        # MCMM_MM3 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_MCMM_MM3.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_MCMM_MM3-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" MCMM    2000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
                f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
                f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_MCMM_MM3".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="SPMC_MM3":
        # SPMC_MM3 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_SPMC_MM3.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_SPMC_MM3-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" SPMC    2000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
                f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
                f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_SPMC_MM3".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="SPMC_OPLS4":
        # SPMC_OPLS4 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_SPMC_OPLS4.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_SPMC_OPLS4-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" SPMC    2000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
                f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
                f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_SPMC_OPLS4".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="LMCS_OPLS4":
        # LMCS_OPLS4 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_LMCS_OPLS4.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_LMCS_OPLS4-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" LMCS    2000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
                f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
                f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
                f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_LMCS_OPLS4".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="LMCS_MM3":
        # LMCS_MM3 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_LMCS_MM3.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_LMCS_MM3-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" LMCS    2000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
                f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
                f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
                f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_LMCS_MM3".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="MCLMCS_MM3":
        # MCLMCS_MM3 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_MCLMCS_MM3.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_MCLMCS_MM3-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" LMCS    2000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
                f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCOP       1      0      0      0     0.5000     0.0000     0.0000     0.0000\n")
                f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
                f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
                f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_MCLMCS_MM3".format(script_directory, name[num])
            os.system(command)
            num+=1
    elif method_name=="MCLMCS_OPLS4":
        # MCLMCS_OPLS4 template
        num=0
        while num < total_number:
            with open(script_directory+name[num]+"_MCLMCS_OPLS4.com", "w+") as f:
                f.write(script_directory+name[num]+"_maestro.mae")
                f.write("\n")
                f.write(script_directory+name[num]+"_MCLMCS_OPLS4-out.mae")
                f.write("\n")
                f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
                f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
                f.write(" LMCS    2000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
                f.write(" NANT       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCNV       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCSS       2      0      0      0    25.0000     0.0000     0.0000     0.0000\n")
                f.write(" MCOP       1      0      0      0     0.5000     0.0000     0.0000     0.0000\n")
                f.write(" DEMX       0    166      0      0    25.0000    50.0000     0.0000     0.0000\n")
                f.write(" COMP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CHIG       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" AUOP       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
                f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
                f.write(" AUTO       0      3      1      1     1.0000     0.0000     0.0000     0.0000\n")
                f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.close()
            command="bmin -WAIT {0}{1}_MCLMCS_OPLS4".format(script_directory, name[num])
            os.system(command)
            num+=1


def CSoutput_split2singlefile():
    command="bash {}results_split.sh".format(script_directory)
    os.system(command)

def frequency_test():
    num=0
    while num < file_count:
        with open(script_directory+minifiles_name[num]+"_fretest.com", "w+") as f:
            f.write(script_directory+minifiles_name[num]+".mae")
            f.write("\n")
            f.write(script_directory+minifiles_name[num]+"_fretest-out.mae")
            f.write("\n")
            f.write(" DEBG     211      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0    500      0     0.0000     0.0000     0.1000     0.0000\n")
            f.write(" ELST       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MTST       0      0      0      0     2.0000     0.0000     0.0000     0.0000\n")
            f.write(" RRHO       0      0      1      0   300.0000     1.0000     2.0000     0.0000\n")
        f.close()
        command="bmin -WAIT {0}{1}_fretest".format(script_directory, minifiles_name[num])
        os.system(command)
        num+=1



script_path=os.path.join(os.getcwd(), sys.argv[1])
script_directory=script_path+"/"
dir_path=script_directory
method=sys.argv[1]
generate_comfile(method)
CSoutput_split2singlefile()


# The pattern to match files like "eX_SPMC_MM3-out-mini-Y.mae" but exclude files with "fretest"
pattern = r"^(?!.*fretest)e\d+.*-out-mini-\d+.mae$"
file_outputmae = list_files_in_directory(dir_path, pattern)
file_count = len(file_outputmae)
print("file_count:", file_count)
print("file_outputmae:", file_outputmae)


minifiles_name=[]
for file in file_outputmae:
    minifiles_name.append(file.split(".")[0])
minifiles_name.sort()
print(minifiles_name)

frequency_test()

