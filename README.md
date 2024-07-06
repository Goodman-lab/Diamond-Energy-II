# Diamond Energy II

Mengman Wei, Jonathan M. Goodman

Diamond Energy II is a systematic conformational searching algorithm operating on a diamond lattice framework, version II.

## Contents
1. Release Notes
2. Requirements and Setup
3. Usage

## Release Notes

Diamond Energy II is an updated version of the algorithm focusing on saccharide-like systems, developed in Python 3.6+ and the RDKit environment.

The original algorithm can be referred to at: [Goodman Lab Diamond Energy](https://github.com/Goodman-lab/Diamond_energy).

## Requirements and Setup

The script is compatible with Python 3.6 or higher and RDKit 2020.09.1. Development and testing were conducted on Linux and MacOS operating systems.

## Getting Started

Ensure Python 3.6+ is installed on your machine and that RDKit is installed and properly configured. If not, installation instructions can be found on the [RDKit website](https://www.rdkit.org/docs/Install.html).

1. Download the `DiamondEnergyII.py` script.
2. Navigate to the directory where the script is located using the terminal.
3. Activate the RDKit environment and run the program from the terminal.

## Correct Usage Syntax

### Global Search Mode

To perform a comprehensive conformational search, use the command:

```bash
python DiamondEnergyII.py '<Molecular_InChI>'

For example, to perform a conformational search for Methylcyclohexane, use the following command:

python DiamondEnergyII.py "InChI=1S/C7H14/c1-7-5-3-2-4-6-7/h7H,2-6H2,1H3"

## Output

The program automatically conducts a conformational search for the specified molecule, such as Methylcyclohexane in the provided example. The terminal will display various details about the search process and its results. The output may include, but is not limited to:

Lowest Energy Conformation: The conformation with the minimum energy.
Energies of Each Conformation: Energy values associated with each identified conformation.
Total Number of Accessible Conformations: Count of conformations that have accessible energy levels.
...and more.
```


## Reference
```bibtex
@phdthesis{Wei2024DiamondEnergy,
  title={*Diamond Energy* – a systematic conformation searching method},
  author={Wei, Mengman},
  year={2024},
  school={University of Cambridge}
  doi={https://doi.org/10.17863/CAM.109200}
}
