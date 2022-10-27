import sys
import gzip

'''The following script converts a MARCS model atmosphere (in a .mod.gz format) into a NICOLE atmosphere file (in a .model format). The script only uses the following columns from the MARCS file: log(tau5000), temperature (in K), electron pressure (in dyn/cm2), microturbulence (in cm/s). The magnetic fields and line-of-sight velocities are assumed to be 0 in the NICOLE file.''' 

'''To run the following script, run python MARCS_to_NICOLE.py <Input file> <Output file>, where the <Input file> is the MARCS model atmosphere (in a .mod.gz format) and <Output file> is the name of the resulting NICOLE atmosphere file'''

MARCS_file = gzip.open(sys.argv[1], "r") #MARCS model file

NICOLE_file = open(sys.argv[2], "w+") #NICOLE model file

NICOLE_file.write("Format version: 1.0\n")
NICOLE_file.write("  1e5    0.\n")

Lines = MARCS_file.readlines()
count = 0
found = False
for line in Lines:
    words = line.split()
    for i in range(0, len(words)):
        words[i] = words[i].decode("utf-8")
    if (words[0] == 'k' and words[1] == 'lgTauR' and words[2] == 'KappaRoss' and found):
        found = False
    if (found):
        for i in range(0, len(words)):
            words[i] = words[i].replace("'", "")
        NICOLE_file.write("  " + str(words[1]) + "  " + str(words[4]) + "  " + str(words[5]) + "  " + str(words[8]) + "  0.000  0.000  0.000   0.000\n")
    if (words[0] == 'k' and words[1] == 'lgTauR' and words[2] == 'lgTau5' and not found):
        found = True
        
MARCS_file.close()
NICOLE_file.close()


