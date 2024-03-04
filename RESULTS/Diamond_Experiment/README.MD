# Diamond Energy II Algorithm Experiments Batch Mode Execution

This repository contains the scripts and files utilized for executing the Diamond Energy II algorithm experiments in batch mode.

## Contents

1. **Test_list.txt**
   - Contains all the test molecules' InChI used in these experiments.

2. **Energyvalue.txt**
   - Introduces the new energy term that is featured in this version of the algorithm.

3. **Sdftomae.sh**
   - Converts file formats during the process. It transforms the conformer output (`.sdf`) from the Diamond Energy II algorithm into a format (``.mae``) that can be used in the MacroModel suite for minimization.

4. **Split_directory.py**
   - The initial script to be run; it generates subdirectories for each test molecule. Each subdirectory is named after the respective test molecule's name. Required files are also automatically copied into each subdirectory. 
   - Conformer search processes are then automatically run for each test molecule within their respective subdirectory.

5. **Generate_mol.py**
   - This script invokes the `diamondenergyII.py` script for each test molecule within their respective subdirectory and manages the conformer search process within that directory.

6. **diamondenergyII.py**
   - Contains the Diamond Energy II algorithm.

7. **Read_energyvalue.py**
   - After the processes are complete, this script reads the energy value results from all subdirectories and compiles them into a single `.csv` file in the main directory.

8. **Read_runningtime.py**
   - Similar to `Read_energyvalue.py`, after completion, it reads the running time results from all subdirectories and compiles them into a single `.csv` file in the main directory.

