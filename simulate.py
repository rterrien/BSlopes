import sys
import subprocess
from datetime import datetime

'''The following script takes inputs from an INPUTS directory, which should contain a text file for the temperature inputs, a text file for the magnetic field inputs, and a text file of the VALD input. The script creates converts MARCS model atmospheres to NICOLE model atmospheres, applies magnetic fields based on the inputs, converts VALD inputs to a NICOLE LINES file, simulaes all the lines using NICOLE, fits the structural parameters of each line using a linear fit, and creates a final table containing relevant information for each line (both VALD inputs and sensitivities to temperature and magnetic field). To run the following script, run 'python simulate.py', assuming the INPUTS directory has the correct inputs. All output information is available in a folder called OUTPUTS/<time>, where <time> is the time the command was executed (to prevent overwriting files)'''

T_input = "INPUTS/temperature_input.txt"
B_input = "INPUTS/magnetic_field_input.txt"
VALD_input = "INPUTS/VALD_input.txt"

subprocess.call(['rm', '-r', 'OUTPUTS'])
subprocess.call(['mkdir', 'OUTPUTS'])
subprocess.call(['cp', '-R', 'INPUTS', 'OUTPUTS'])

now = datetime.now()
year = now.year
month = now.month
day = now.day
hour = now.hour
minute = now.minute
sec = now.second
time = str(month) + "-" + str(day) + "-" + str(year) + "-" + str(hour) + ":" + str(minute) + ":" + str(sec)

output = 'OUTPUT_' + str(time)

subprocess.call(['python', 'create_model_atmospheres.py', T_input, B_input])

subprocess.call(['python', 'VALD_to_NICOLE.py', VALD_input])

subprocess.call(['python', 'run_VALD_lines.py', T_input, B_input, 'OUTPUTS/LINES'])

subprocess.call(['python', 'line_fitting.py', T_input, B_input, 'OUTPUTS/LINES'])

subprocess.call(['python', 'create_table.py', 'OUTPUTS/SLOPES.txt', 'OUTPUTS/VALD_TABLE.txt'])

subprocess.call(['mv', '-v', 'OUTPUTS', output])






