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


# The molecular name in extension work test set
# One can change the following list to target test set
name=["e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9",
      "e10", "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19",
      "e20", "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29",
      "e30", "e31", "e32", "e33", "e34", "e35", "e36", "e37"]
total_number=len(name)


# The comparing experiment names of Macromodel with algorithm and force filed
Maestro_test=["MCMM_OPLS4", "MCMM_MM3", "SPMC_OPLS4", "SPMC_MM3", "LMCS_OPLS4", "LMCS_MM3", "MCLMCS_OPLS4", "MCLMCS_MM3"]


#def count_files_in_directory(directory, extension):
#    return len([f for f in os.listdir(directory) if f.endswith(extension)])
#dir_path=os.getcwd()
#file_extension='.sdf'
#file_count=count_files_in_directory(dir_path, file_extension)
#print(file_count)
def count_files_in_directory(directory, extension):
    return len([f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith(extension.lower())])

def files_in_directory(directory, extension):
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith(extension.lower())]




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
            f.write(" MINI       1      0  10000      0     0.0000     0.0000     0.0000     0.0000\n")
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
            f.write(" MINI       1      0  10000      0     0.0000     0.0000     0.0000     0.0000\n")
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
            f.write(" MINI       1      0  10000      0     0.0000     0.0000     0.0000     0.0000\n")
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
            f.write(" MINI       1      0  10000      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        num+=1
'''


# Submit batch MacroModel CS running
def run_CS_batch():
    a=0
    while a < total_number:
        command="bmin "+name[a]+"_MCMM_OPLS4"
        os.system(command)
        a+=1


def run_diamond():
    # need to adjust for batch run later on
    #command="python "+dir_path+"diamond_oxy8.py '"+molinchi+"' "+"0 "+molname
    #os.system(command)
    # Change to the desired directory
    os.chdir(dir_path)

    # Execute the command
    command = "python diamond_oxy8.py '" + molinchi + "' 0 " + molname
    os.system(command)


def Add_hydrogen():
    c=0
    while c < file_count:
        suppl_mol=Chem.SDMolSupplier(script_directory+file_addHs[c]+".sdf")
        mol=suppl_mol[0]
        mol_addhs=Chem.AddHs(mol, addCoords=True)
        # Do not use embed as it will change the carbon chain's structure losing diamond output feature
        #AllChem.EmbedMolecule(mol_addhs)
        mol_addhs_file=Chem.MolToMolFile(mol_addhs, script_directory+file_addHs[c]+"H.sdf")
        c+=1


def call_sdftomae():
    command="bash {}sdftomae.sh".format(script_directory)
    os.system(command)

import glob

def mae_addHs():
    for mae_file in glob.glob("*.mae"):
        # Skip files that have already been treated (i.e., filenames ending with 'H.mae')
        if mae_file.endswith("H.mae"):
            continue

        output_file = mae_file.replace('.mae', 'H.mae')
        command = '/shared/shared/schrodinger/2022-1/utilities/applyhtreat {} {} -t "All-atom with No-Lp"'.format(mae_file, output_file)
        os.system(command)


# current OPLS4, change to MM3 for the same used
def mini4diamond():
    num=0
    while num < file_count:
        with open(script_directory+minifiles_name[num]+"_minimization.com", "w+") as f:
            f.write(script_directory+minifiles_name[num]+".mae")
            f.write("\n")
            f.write(script_directory+minifiles_name[num]+"_minimization-out.mae")
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
            f.write(" MINI       1      0   2000      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
        f.close()
        command = "bmin" + " -WAIT " + os.path.join(script_directory, minifiles_name[num] + "_minimization")
        os.system(command)
        num+=1

def frequency_test():
    num=0
    while num < file_count:
        with open(script_directory+minifiles_name[num]+"_fretest.com", "w+") as f:
            f.write(script_directory+minifiles_name[num]+"_minimization-out.mae")
            f.write("\n")
            f.write(script_directory+minifiles_name[num]+"_fretest-out.mae")
            f.write("\n")
            f.write(" DEBG     211      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n")
            f.write(" BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n")
            f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MINI       4      0    500      0     0.0000     0.0000     0.1000     0.0000\n")
            f.write(" ELST       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
            f.write(" MTST       0      0      0      0     2.0000     0.0000     0.0000     0.0000\n")
            f.write(" RRHO       0      0      1      0   300.0000     1.0000     2.0000     0.0000\n")
        f.close()
        command = "bmin" + " -WAIT " + os.path.join(script_directory, minifiles_name[num] + "_fretest")
        os.system(command)
        num+=1



def list_files_in_directory(directory, extension, exclude_substring=None):
    if exclude_substring:
        return [f for f in os.listdir(directory) if f.endswith(extension) and exclude_substring not in f]
    else:
        return [f for f in os.listdir(directory) if f.endswith(extension)]


#inchi_to_sdf()
#generate_comfile()
#run_CS_batch()
script_path=os.path.join(os.getcwd(), sys.argv[1])
script_directory=script_path+"/"
dir_path=script_directory
#script_directory=os.getcwd()+"/"+sys.argv[1]+"/"
#print("script_directory:", script_directory)
#molinchifile=[f for f in script_directory if f.endswith('.txt')]
molinchifile = list_files_in_directory(dir_path, ".txt", "data") 
print(molinchifile)
print("molinchifile:", molinchifile)
print("length:", len(molinchifile))
print("type:", type(molinchifile))
molname=molinchifile[0].split(".")[0]
print("molname:", molname)
with open(script_directory+molinchifile[0], "r") as f:
    molinchi=f.readline().split("\n")[0]
f.close()


run_diamond()
file_extension='.sdf'
file_count=count_files_in_directory(dir_path, file_extension)
print("dir_path:", dir_path)
print("file_count:", file_count)

file_outputsdf=files_in_directory(dir_path, file_extension)
#file_addHs=[]
#for file in file_outputsdf:
#    file_addHs.append(file.split(".")[0])
#file_addHs.sort()
#print(file_addHs)
#Add_hydrogen()

call_sdftomae()
mae_addHs()

# count and collect the input mol files end with .mae for further minimization and frequency test
file_outputmae=list_files_in_directory(dir_path, "H.mae", "-out") # file pattern
minifiles_name=[]
for file in file_outputmae:
    minifiles_name.append(file.split(".")[0])

mini4diamond()
frequency_test()
