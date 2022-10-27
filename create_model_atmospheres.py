import sys
import subprocess

'''The following script creates a list of NICOLE model atmosphere files based on a list of temperatures and magnetic fields (given as textfiles)'''

'''To run the following script, run "python create_model_atmospheres.py <Temperature input file> <Magnetic Field input file>"where <Temperature input file> is a text file, with each line containing a temperature value and a MARCS model atmosphere (separated with a space) and <Magnetic Field input file> is a text file, with each line containing a magnetic field value'''
T_input = open(sys.argv[1], "r")
T_lines = T_input.readlines()
B_input = open(sys.argv[2], "r")
B_lines = B_input.readlines()
T_list = []
B_list = []
MARCS_model_atmosphere_list = []

#Create a directory that will contain all NICOLE model atmosphere files
subprocess.call(['mkdir', 'OUTPUTS/model_atmospheres'])

#Read tempature values and MARCS model atmospheres from the temperature input file
for i in T_lines:
    words = i.split()
    T_list.append(words[0])
    MARCS_model_atmosphere_list.append("INPUTS/" + str(words[1]))

#Read magnetic field values from the magnetic field input file
for j in B_lines:
    words = j.split()
    B_list.append(words[0])

#Convert each MARCS model atmosphere to a NICOLE model atmosphere
for k in range(0, len(T_list)):
    subprocess.call(['python', 'MARCS_to_NICOLE.py', MARCS_model_atmosphere_list[k], "OUTPUTS/model_atmospheres/T" + str(T_list[k]) + "_MF0.model"])

#Create all NICOLE model atmospheres for all temperatures and magnetic fields
for i in T_list:
    for j in B_list:
        subprocess.call(['python', 'convert.py', 'OUTPUTS/model_atmospheres/T' + str(i) + '_MF0.model', 'OUTPUTS/model_atmospheres/T' + str(i) + '_MF' + str(j) + '.model', 'x=' + str(j), 'y=' + str(j), 'z=' + str(j)])           
