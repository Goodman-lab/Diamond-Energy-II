import os
import re
import pandas as pd


def convert_time_to_seconds(time_str):
    # Remove the trailing 's' if it exists
    time_str = time_str.rstrip('s')

    if 'm' in time_str:
        mins, secs = time_str.split('m')
        return float(mins) * 60 + float(secs)
    else:
        return float(time_str)


# List of molecular names and test methods
name=["e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9",
      "e10", "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19",
      "e20", "e21", "e22", "e23", "e24", "e25", "e26", "e27", "e28", "e29",
      "e30", "e31", "e32", "e33", "e34", "e35", "e36", "e37"]

Maestro_test = ["MCMM_OPLS4", "MCMM_MM3", "SPMC_OPLS4", "SPMC_MM3",
                "LMCS_OPLS4", "LMCS_MM3", "MCLMCS_OPLS4", "MCLMCS_MM3"]

def read_eachMaestrotest():
    for test_method in Maestro_test:
        final_df = pd.DataFrame(columns=['Molecule', 'Unique_Structures', 'Generated_Structures', 'Minimization_Time', 'Total_running_time', 'conformer_generate_time', 'ave_conformer_generate_time', 'ave_running_time'])

        for molecule in name:
            file_path = os.path.join(".", test_method, f"{molecule}_{test_method}.log")
            try:
                with open(file_path, "r", encoding='utf-8') as f:
                    lines = f.readlines()

                    Output_number = None
                    generated_structures = None
                    conformer_minimization_time = None
                    total_running_time = 0
                    ave_conformer_generated_time = 0
                    ave_running_time = 0

                    for line in lines:
                        line = line.strip()

                        if "minimized with good convergence" in line:
                            try:
                                Output_number = int(line.split("minimized")[0].strip())
                            except ValueError:
                                continue

                        generated_match = re.search(r'(\d+)\s+structures generated', line)
                        if generated_match:
                            generated_structures = int(generated_match.group(1))

                        time_match = re.search(r'Time in energy minimizations:\s+([\d.]+)\s+CPU sec', line)
                        if time_match:
                            conformer_minimization_time = float(time_match.group(1))

                        # Capture user time in both formats
                        # Capture total running time (user's CPU time)
                        total_time_match = re.search(r'([\d.]+)user', line)
                        if total_time_match:
                            total_running_time += float(total_time_match.group(1))
                            continue  # Skip to the next line


                        total_time_match = re.search(r'user\s*([\d\.:m]+[s]*)', line)
                        if total_time_match:
                            total_running_time = convert_time_to_seconds(total_time_match.group(1))
                            continue  # Skip to the next line

                    if total_running_time is not None and conformer_minimization_time is not None:
                        conformer_generate_time = total_running_time - conformer_minimization_time
                        ave_conformer_generated_time = conformer_generate_time / generated_structures
                        ave_running_time = total_running_time / Output_number
                    else:
                        conformer_generate_time = None

                    new_row = pd.DataFrame({
                        'Molecule': [molecule],
                        'Unique_Structures': [Output_number],
                        'Generated_Structures': [generated_structures],
                        'Minimization_Time': [conformer_minimization_time],
                        'Total_running_time': [total_running_time],
                        'conformer_generate_time': [conformer_generate_time],
                        'ave_conformer_generate_time': [ave_conformer_generated_time],
                        'ave_running_time': [ave_running_time]
                    })

                    final_df = pd.concat([final_df, new_row], ignore_index=True)

            except FileNotFoundError:
                print(f"File {file_path} not found.")
                continue

        final_df.to_csv(os.path.join(".", test_method, '_total_running_time.csv'), index=False)

# Execute the function to read the log files
read_eachMaestrotest()
