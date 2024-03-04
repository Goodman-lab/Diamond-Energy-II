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

# The molecular name in extension work test set
# One can change the following list to target test set
name=["e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9",
      "e10", "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19",
      "e20", "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29",
      "e30", "e31", "e32", "e33", "e34", "e35", "e36", "e37"]
total_number=len(name)


# The complete test molecule InChi information
with open("Test_list.txt", "r") as f:
    inchi_list=[]
    for line in f.read().split():
        print(line)
        inchi_list.append(line)


def count_files_in_directory(directory, extension):
    return len([f for f in os.listdir(directory) if f.endswith(extension)])

#dir_path='.'
#file_extension='.mae'
#file_count=count_files_in_directory(dir_path, file_extension)

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
            AllChem.EmbedMolecule(mol_maestro)
            mol_maetsro_file=Chem.MolToMolFile(mol_maestro, name[inchi_n]+"_maestro.sdf")
            inchi_n+=1

# Generate input mol file for running MacroModel
# The function can use structconvert.py in Schrodinger/version/utilities/ directory
# To change .sdf format mol file to .mae MacroModel molecular structure format as input
# More instructions can be found in macromodel_reference_manual
def sdf2mae():
    i=0
    while i < total_number:
        # The below directory to call structconvert.py may vary in different computer
        command="/shared/shared/schrodinger/2022-1/utilities/structconvert "+name[i]+".sdf "+name[i]+"_maestro.mae"
        os.system(command)
        i+=1

# Generate command file for running MacroModel
# to run MacroModel from command line for batch running, instead of using its graphic interaction software Maestro by hand
# We need a molecular structure file in .mae format, which can be converted using structconvert.py in Schrodinger/version/utilities/ directory
# and a command file in .com format, which includes all the tasks (operation codes and instructions can be found in macromodel_reference_manual)
# function for generating MacroModel command .com file for corresponding molecular structure .mae file

# MCMM_OPLS4 template
def generate_comfile():
    num=0
    while num < total_number:
        with open(name[num]+"_MCMM_OPLS4.com", "w+") as f:
            f.write(name[num]+"_maestro.mae")
            f.write("\n")
            f.write(name[num]+"_MCMM_OPLS4-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" MCMM  100000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
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
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1
'''
# SPMC_MM3 template
def generate_comfile():
    num=0
    while num < total_number:
        with open(name[num]+"_SPMC_MM3.com", "w+") as f:
            f.write(name[num]+"_maestro.mae")
            f.write("\n")
            f.write(name[num]+"_SPMC_MM3-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" SPMC  100000      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
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
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1



# LMCS_OPLS4 template
def generate_comfile():
    num=0
    while num < total_number:
        with open(name[num]+"_LMCS_OPLS4.com", "w+") as f:
            f.write(name[num]+"_maestro.mae")
            f.write("\n")
            f.write(name[num]+"_LMCS_OPLS4-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" LMCS  100000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
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
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1
        
# MCLMCS_MM3 template
def generate_comfile():
    num=0
    while num < total_number:
        with open(name[num]+"_MCLMCS_MM3.com", "w+") as f:
            f.write(name[num]+"_maestro.mae")
            f.write("\n")
            f.write(name[num]+"_MCLMCS_MM3-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" LMCS  100000      0      0      0     0.0000     0.0000     3.0000     6.0000\n")
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
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1
'''

# Submit batch MacroModel CS running
def run_CS_batch():
    a=0
    while a < total_number:
        command="bmin -WAIT "+name[a]+"_MCMM_OPLS4"
        os.system(command)
        a+=1

def run_diamond():
    # need to adjust for batch run later on
    command="python Diamond_energy.py 'InChI=1S/C18H38/c1-3-5-7-9-11-13-15-17-18-16-14-12-10-8-6-4-2/h3-18H2,1-2H3'"
    os.system(command)

def call_sdftomae():
    command="bash sdftomae.sh"
    os.system(command)

def generate_minifile4diamond():
    num=0
    while num < file_count:
        with open(str(num)+"_diamond_minimization.com", "w+") as f:
            f.write(str(num)+"_diamond.mae")
            f.write("\n")
            f.write(str(num)+"_diamond_minimization-out.maegz")
            f.write("\n")
            f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG      55      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" DEBG    1003      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" CRMS       0      0      0      0     0.0000     0.5000     0.0000     2.0000\n")
            f.write(" BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" CONV       2      0      0      0     0.0010     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0   1000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1

# Submit batch minimization running for Diamond_output
def run_mini_batch():
    b=0
    while b < total_number:
        command="bmin "+str(b)+"_diamond_minimization"
        os.system(command)
        b+=1

def split_directory():
    mol_order=0
    while mol_order < total_number:
        command="mkdir "+name[mol_order]
        os.system(command)
        #Best write it down and then can check the structure later on
        with open("./"+name[mol_order]+"/"+name[mol_order]+".txt", "w+") as f:
            f.write(inchi_list[mol_order])
            f.write("\n")
        f.close()
        command1="cp diamondenergyII.py ./"+name[mol_order]
        os.system(command1)
        command2="cp generate_mol.py ./"+name[mol_order]
        os.system(command2)
        command3="cp sdftomae.sh ./"+name[mol_order]
        os.system(command3)
        command4="cp energyvalue.txt ./"+name[mol_order]
        os.system(command4)
        mol_order+=1


def run_diamondscript_parallel():
    mol_test=0
    while mol_test < total_number:
        current_running_path="./"+name[mol_test]+"/"
        command_run="nohup python -u "+current_running_path+"generate_mol.py "+name[mol_test]+"> "+current_running_path+"running_record.log 2>&1 &"
        os.system(command_run)
        mol_test+=1



#inchi_to_sdf()
#generate_comfile()
#run_CS_batch()
#run_diamond()
#call_sdftomae()
#generate_minifile4diamond()
#run_mini_batch()
#split_directory()
run_diamondscript_parallel()
