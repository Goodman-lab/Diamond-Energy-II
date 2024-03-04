# MacroModel Conformer Search Algorithms Batch Execution

This repository contains scripts and files for running batch mode experiments with four leading conformer search algorithms (SPMC, LMCS, MCLMCS, MCMM) in the MacroModel suite. These scripts require a valid Schrödinger license to execute the models within the MacroModel suite.

## Contents

1. **Test_list.txt**
   - Lists all the test molecules' InChI codes used in these experiments.

2. **Results_split.sh**
   - Splits the compact minimization results (`.maegz`) for a given test molecule into individual files (`.mae`), which can then be used for frequency testing in MacroModel.

3. **Split_directory_mae.py**
   - The initial script to be run; it creates subdirectories for each test algorithm, naming each one after the algorithm name and the force field used. Required files are also automatically copied into each subdirectory.
   - The conformer search process is then automatically executed for all test molecules in their corresponding subdirectory, using the specified conformer search method and force field.

4. **maestro_mol.py**
   - Executes each conformer search method and force field test for all test molecules within their respective subdirectory and manages the conformer search process within that directory.

5. **Read_energyvalue.py**
   - Once all processes are complete, this script is used to compile energy value results from all subdirectories into a single `.csv` file in the main directory.

6. **Read_runningtime.py**
   - After completion, this script compiles running time results from all subdirectories into a single `.csv` file in the main directory.
