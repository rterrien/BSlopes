import sys
from datetime import datetime
import subprocess

'''The following script runs all lines in the LINES file using NICOLE. It outputs all .pro files in a directory named by the current time of the simulations.'''

'''To run the following script, run "python run_VALD_lines.py <Temperature input file> <Magnetic Field input file> <LINES>" where <Temperature input file> is a text file, with each line containing a temperature value and a MARCS model atmosphere (separated with a space), <Magnetic Field input file> is a text file, with each line containing a magnetic field value, and <LINES> is a textfile containing a list of LINES. The script assumes all model atmosphere files are in a directoy called 'OUTPUTS/model_atmospheres' '''

T_input = open(sys.argv[1], "r")
T_lines = T_input.readlines()
B_input = open(sys.argv[2], "r")
B_lines = B_input.readlines()
LINES_input = sys.argv[3]
subprocess.call(['cp', LINES_input, 'LINES'])

T_list = []
B_list = []
MARCS_model_atmosphere_list = []

#Read tempature values and MARCS model atmospheres from the temperature input file
for i in T_lines:
    words = i.split()
    T_list.append(words[0])
    MARCS_model_atmosphere_list.append(words[1])

#Read magnetic field values from the magnetic field input file
for j in B_lines:
    words = j.split()
    B_list.append(words[0])

#Read the LINES file and extract the species and wavelength of each line
VALD_LINES_file = open(sys.argv[3], "r")
VALD_LINES = VALD_LINES_file.readlines()
species = []
wavelengths = []
for line in VALD_LINES:
    words = line.split()
    if len(words) > 0 and words[0][0] == '[' and words[-1][-1] == ']':
        species.append(words[0][1::])
        wavelengths.append(words[-1][0:-1])

temperatures = T_list
magnetic_fields = B_list

#Make directories for all Figures
subprocess.call(['mkdir', 'OUTPUTS/model_files'])
subprocess.call(['mkdir', 'OUTPUTS/Figures'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/varying_B'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/varying_T'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/EWvsT'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/EWvsB'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/FWHMvsT'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/FWHMvsB'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/DepthvsT'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/DepthvsB'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/EWfvsT'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/EWfvsB'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/FWHMfvsT'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/FWHMfvsB'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/DepthfvsT'])
subprocess.call(['mkdir', 'OUTPUTS/Figures/DepthfvsB'])


#Open the NICOLE.input file
NICOLE_file = open("NICOLE.input", "r")
NICOLE_lines = NICOLE_file.readlines()
new_NICOLE_lines = []

#Iterate over each line
for k in range(0, len(species)):
    element = species[k]
    central_wavelength = wavelengths[k]
    print(central_wavelength)
    
    #Iterate over each temperature
    for i in temperatures:

        #The control model file for varying temperature is defined as the atmosphere with the lowest (0 G) magnetic field
        model_file = "OUTPUTS/model_atmospheres/T" + str(i) + "_MF0.model"
        
        #Read the NICOLE.input file and read the model atmosphere file, the wavelength paramaters, and the LINE name
        for line in NICOLE_lines:
            words = line.split()
            if (len(words) >= 3 and words[0] == "Input" and words[1] == "model="):
                input_model_name = model_file
                new_NICOLE_lines.append("Input model= " + str(model_file) + "\n")
            elif (len(words) >= 3 and words[0] == "First" and words[1] == "wavelength="):
                first_lambda = str(float(central_wavelength) - 2)
                new_NICOLE_lines.append("  First wavelength= " + str(first_lambda) + "\n")
            elif (len(words) >= 3 and words[0] == "Wavelength" and words[1] == "step="):
                lambda_step = words[2]
                new_NICOLE_lines.append(line)
            elif (len(words) >= 4 and words[0] == "Number" and words[1] == "of" and words[2] == "wavelengths="):
                num_lambda = words[3]
                new_NICOLE_lines.append(line)
            elif (len(words) >= 2 and words[0][0:5] == "Line="):
                new_NICOLE_lines.append("Line=" + str(element) + " " + str(central_wavelength) + "\n")
            else:
                new_NICOLE_lines.append(line)
        file_name = str(element) + "_" + str(float(central_wavelength)) + "_T" + str(i) + "_" + "MF0" + "_" + lambda_step + "_" + num_lambda
 
        #Read the NICOLE.input file and write the new names of the .mod and .pro files
        output = open("NICOLE.input", "w+")
        for line in new_NICOLE_lines:
            words = line.split()
            if (len(words) >= 3 and words[0] == "Output" and words[1] == "profiles="):
                new_line = "Output profiles= OUTPUTS/model_files/" + str(file_name) + ".pro # Output profiles \n"
                output.write(new_line)
            elif (len(words) >= 3 and words[0] == "Output" and words[1] == "model="):
                new_line = "Output model= OUTPUTS/model_files/" + str(file_name) + ".mod # Output profiles \n"
                output.write(new_line)
            else:
                output.write(line)
        output.close()
        subprocess.call(['python', 'run_nicole.py'])
        new_NICOLE_lines = []


    #Iterate over each temperature
    for i in magnetic_fields:

        #The control model file for varying magnetic field is defined as the atmosphere with the lowest temperature
        model_file = "OUTPUTS/model_atmospheres/T" + str(temperatures[0]) + "_MF" + str(i) + ".model"
        
        #Read the NICOLE.input file and read the model atmosphere file, the wavelength paramaters, and the LINE name
        for line in NICOLE_lines:
            words = line.split()
            if (len(words) >= 3 and words[0] == "Input" and words[1] == "model="):
                input_model_name = model_file
                new_NICOLE_lines.append("Input model= " + str(model_file) + "\n")
            elif (len(words) >= 3 and words[0] == "First" and words[1] == "wavelength="):
                first_lambda = str(float(central_wavelength) - 2)
                new_NICOLE_lines.append("  First wavelength= " + str(first_lambda) + "\n")
            elif (len(words) >= 3 and words[0] == "Wavelength" and words[1] == "step="):
                lambda_step = words[2]
                new_NICOLE_lines.append(line)
            elif (len(words) >= 4 and words[0] == "Number" and words[1] == "of" and words[2] == "wavelengths="):
                num_lambda = words[3]
                new_NICOLE_lines.append(line)
            elif (len(words) >= 2 and words[0][0:5] == "Line="):
                new_NICOLE_lines.append("Line=" + str(element) + " " + str(central_wavelength) + "\n")
            else:
                new_NICOLE_lines.append(line)
        file_name = str(element) + "_" + str(float(central_wavelength)) + "_T" + str(temperatures[0]) + "_MF" + str(i) +  "_" + lambda_step + "_" + num_lambda

        #Read the NICOLE.input file and write the new names of the .mod and .pro files
        output = open("NICOLE.input", "w+")
        for line in new_NICOLE_lines:
            words = line.split()
            if (len(words) >= 3 and words[0] == "Output" and words[1] == "profiles="):
                new_line = "Output profiles= OUTPUTS/model_files/" + str(file_name) + ".pro # Output profiles \n"
                output.write(new_line)
            elif (len(words) >= 3 and words[0] == "Output" and words[1] == "model="):
                new_line = "Output model= OUTPUTS/model_files/" + str(file_name) + ".mod # Output profiles \n"
                output.write(new_line)
            else:
                output.write(line)
        output.close()
        subprocess.call(['python', 'run_nicole.py'])
        new_NICOLE_lines = []

