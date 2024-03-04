import os
import re
import pandas as pd

def read_running_time():
    final_df = pd.DataFrame(columns=['Molecule', 'Geometry_generate_time', 'Minimization_time', 'Total_Running_Time', 'Average_Running_Time', 'Total_conformer_num', 'Output_conformer_num'])

    # Navigate through each subdirectory
    for dirpath, dirnames, filenames in os.walk('.'):
        geometry_time = 0
        total_time = 0
        minimization_count = 0
        total_conformer_num = 0
        output_conformer_num = 0

        # Check each file in the current subdirectory
        for filename in filenames:
            if filename.endswith("minimization.log"):
                minimization_count += 1
                with open(os.path.join(dirpath, filename), "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        line = line.strip()
                        # Capture total running time (user's CPU time)
                        total_time_match = re.search(r'([\d.]+)user', line)
                        if total_time_match:
                            total_time += float(total_time_match.group(1))

            elif filename == "running_record.log":
                with open(os.path.join(dirpath, filename), "r") as f:
                    lines = f.readlines()
                    for line in lines:
                        if "Diamond_running_time:" in line:
                            geometry_time = float(line.split(":")[1].strip())
                        elif "Number of conformations with accessible energies" in line:
                            output_conformer_num = int(line.split(":")[1].strip())
                        elif "rotatable bonds;" in line:
                            total_conformer_num = int(line.split(";")[1].strip().split()[0])

        total_running_time = total_time + geometry_time
        average_time = total_running_time / output_conformer_num if output_conformer_num != 0 else 0

        # Create a new row as a DataFrame and concatenate it with the final DataFrame
        new_row = pd.DataFrame([{
            'Molecule': os.path.basename(dirpath),
            'Geometry_generate_time': geometry_time,
            'Minimization_time': total_time,
            'Total_Running_Time': total_running_time,
            'Average_Running_Time': average_time,
            'Total_conformer_num': total_conformer_num,
            'Output_conformer_num': output_conformer_num
        }])
        final_df = pd.concat([final_df, new_row], ignore_index=True)
    # Save the DataFrame as a CSV file
    final_df.to_csv('Diamond_OPLS4_MacroModel_total_running_times.csv', index=False, float_format='%.15f')

read_running_time()

