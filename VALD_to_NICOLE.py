import sys
import subprocess
import math

'''The following script converts a list of VALD lines (using the long extraction format) into a NICOLE LINES file. The script also outputs the parameters in a table format, a textfile of the lines that do not conform to the proper format, and a textfile of the duplicate lines'''

'''To run the following script, run "python VALD_to_NICOLE.py <Input file>" where the <Input file> is  the VALD linelist, extracted in the long format'''

VALD_output = open(sys.argv[1], "r")
VALD_lines = VALD_output.readlines()
LINES_file = open('OUTPUTS/LINES', 'w')
error_file = open('OUTPUTS/ERRORS.txt', 'w')
TABLE = open('OUTPUTS/VALD_TABLE.txt', 'w')
duplicates_file = open('OUTPUTS/DUPLICATES.txt', 'w')
first_word = []
element_symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th,' 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']

alphabet = list(map(chr, range(97, 123)))
alphabet = [x.upper() for x in alphabet]

#Function to convert integers to Roman numerals
def intToRoman(num):
  
    # Storing roman values of digits from 0-9
    # when placed at different places
    m = ["", "M", "MM", "MMM"]
    c = ["", "C", "CC", "CCC", "CD", "D",
         "DC", "DCC", "DCCC", "CM "]
    x = ["", "X", "XX", "XXX", "XL", "L",
         "LX", "LXX", "LXXX", "XC"]
    i = ["", "I", "II", "III", "IV", "V",
         "VI", "VII", "VIII", "IX"]
  
    # Converting to roman
    thousands = m[num // 1000]
    hundreds = c[(num % 1000) // 100]
    tens = x[(num % 100) // 10]
    ones = i[num % 10]
  
    ans = (thousands + hundreds +
           tens + ones)
  
    return ans

#Define flags and lists of parameters
data_flag = False
error_flag = False
total_line_count = 0
total_error_count = 0
first_line = VALD_lines[0].split(',')
NUM_LINES = int(first_line[2].strip())
all_species = []
all_wavelengths = []
all_elements = []
all_ionizations = []
all_ionization_stages = []
all_excitation_potentials = []
all_log_gfs = []
all_term_lowers = []
all_term_uppers = []
all_collisions = []
all_gamma_radiatives = []
all_gamma_starks = []
all_gamma_van_der_waals = []
all_widths = []
all_depths = []
all_lande_factors = []
spec_ions = []
wavelengths = []
log_gfs = []
E_lows = []
J_lows = []
E_ups = []
J_ups = []
lower_landes = []
upper_landes = []
mean_landes = []
rads = []
starks = []
waalses = []
depths = []
line2s = []
line3s = []
line4s = []

#Iteratively read each line of the VALD output and extract the relevant parameters of each line
for k in range(3, NUM_LINES*4+3): 

    line = VALD_lines[k]
    words = line.split()

    #Parses every first line and extracts parameters such as species, wavelength, etc.
    if k % 4 == 3 and words[0].replace('\'', '') in element_symbols:

        data_flag = True
        spec_ion = words[0]
        spec_ion = spec_ion.replace('\'', '').replace(',', '')
        ionization_num = words[1].replace('\'', '').replace(',', '')
        element = spec_ion       
 
        #Convert ionization number
        if (ionization_num).isnumeric():
            ionization = intToRoman(int(ionization_num))
        else:
            error_flag = True
            error_file.write(VALD_lines[k])
            error_file.write(VALD_lines[k + 1])
            error_file.write(VALD_lines[k + 2])
            total_error_count += 1

        wavelength = words[2].replace(',', '')
        species_name = element + ionization_num
        log_gf = words[3].replace(',', '')
        E_low = words[4].replace(',', '')
        J_low = words[5].replace(',', '')
        E_up = words[6].replace(',', '')
        J_up = words[7].replace(',', '')
        lower_lande = words[8].replace(',', '')
        upper_lande = words[9].replace(',', '')
        mean_lande = words[10].replace(',', '')
        excitation_potential = words[4].replace(',', '') # in ev
        lower_term_2 = words[5].replace(',', '')
        upper_term_2 = words[7].replace(',', '')
        lande_factor = words[10].replace(',', '')
        depth = words[-1].replace(',', '')
        
        #Parse the gamma parameters
        if len(words) == 13 and len(words[11].split(',')) >= 3:
            collisions = 1 #should be collisions = 3 if using gamma paramaters
            gammas = words[11].split(',')
            rad = gammas[0]
            stark = gammas[1]
            waals = gammas[2]
            #convert the units of gamma parameters to fit NICOLE's requirements
            #rad = round(4 * math.pi * 10**float(rad) / 10**8, 2)
            #stark = round(4 * math.pi * 10**float(stark) * 10**4, 2)
            #waals = round(4 * math.pi * 10**float(waals) * 10**8, 2)
        elif len(words) >= 14:
            collisions = 1
            rad = words[11].replace(',', '')
            stark = words[12].replace(',', '')
            waals = words[13].replace(',', '')
        else:
            collisions = 1
            rad = 0.000
            stark = 0.000
            waals = 0.000
            error_file.write(VALD_lines[k])
            error_file.write(VALD_lines[k + 1])
            error_file.write(VALD_lines[k + 2])
            error_flag = True
            total_error_count += 1

    #Parses every second line and extracts the term designation
    elif k % 4 == 0 and data_flag == True:
        mult_letters = False
        any_letters = any(c.isalpha() for c in words[-1][0:2])
        count = 0
        for j in words[-1]:
            if j in alphabet:
                count += 1
        if count > 1:
            mult_letters = True
        #get rid of badly formatted terms        
        if (not any_letters) or mult_letters or '10' in words[-1] or '=' in words[-1] or '.' in words[-1] or '(' in words[-1] or ')' in words[1] or '<' in words[-1] or '>' in words[-1]:
            error_file.write(VALD_lines[k - 1])
            error_file.write(VALD_lines[k])
            error_file.write(VALD_lines[k + 1])
            error_flag = True
            total_error_count += 1
        else:
            lower_term_1 = ''
            digit_flag = False
            for letter in words[-1]:
                if letter.isdigit():
                    digit_flag = True
                if digit_flag:
                    lower_term_1 += letter
            lower_term_1 = lower_term_1.replace('\'', '').replace('\\', '').replace('(', '').replace(')', '').replace('*', '').upper()
            line2 = line.replace('\'', '')
    
    #Parses every third line and extracts the term designation
    elif k % 4 == 1 and data_flag == True:
        mult_letters = False
        any_letters = any(c.isalpha() for c in words[-1][0:2])
        count = 0
        for j in words[-1]:
            if j in alphabet:
                count += 1
        if count > 1:
            mult_letters = True
        #get rid of badly formatted terms
        if (not any_letters) or mult_letters or '10' in words[-1] or '=' in words[-1] or '.' in words[-1] or '(' in words[-1] or ')' in words[1] or '<' in words[-1] or '>' in words[-1]:
            error_file.write(VALD_lines[k - 2])
            error_file.write(VALD_lines[k - 1])
            error_file.write(VALD_lines[k])
            error_flag = True
            total_error_count += 1
        elif error_flag == False:
            upper_term_1 = ''
            digit_flag = False
            for letter in words[-1]:
                if letter.isdigit():
                    digit_flag = True
                if digit_flag:
                    upper_term_1 += letter
            upper_term_1 = upper_term_1.replace('\'', '').replace('\\', '').replace('(', '').replace(')', '').replace('*', '').upper()
            line3 = line.replace('\'', '')
            
            #append all parameters to their respective lists
            all_species.append(species_name)
            all_wavelengths.append(wavelength)
            all_elements.append(element)
            all_ionizations.append(ionization)
            all_excitation_potentials.append(excitation_potential)
            all_ionization_stages.append(ionization_num)
            all_log_gfs.append(log_gf)
            all_term_lowers.append(lower_term_1 + lower_term_2)
            all_term_uppers.append(upper_term_1 + upper_term_2)
            all_collisions.append(1)
            all_gamma_radiatives.append(rad)
            all_gamma_starks.append(stark)
            all_gamma_van_der_waals.append(waals)
            all_widths.append(2)
            all_depths.append(depth)
            all_lande_factors.append(lande_factor)
            total_line_count += 1
   
    #Reset the flags and ignore the fourth line
    elif data_flag == True and k % 4 == 2:

        if error_flag == False:
            line4 = line.replace('\'', '')

            #append all parameters to their respective lists
            spec_ions.append(spec_ion)
            wavelengths.append(wavelength)
            log_gfs.append(log_gf)
            E_lows.append(E_low)
            J_lows.append(J_low)
            E_ups.append(E_up)
            J_ups.append(J_up)
            lower_landes.append(lower_lande)
            upper_landes.append(upper_lande)
            mean_landes.append(mean_lande)
            rads.append(rad)
            starks.append(stark)
            waalses.append(waals)
            depths.append(depth)
            line2s.append(line2)
            line3s.append(line3)
            line4s.append(line4)

        error_flag = False
        data_flag = False   
        duplicate_flag = False 

#Remove any duplicate lines with the same species and wavelength. Only keep the lines with the larger depth
i = 0
count = 0
duplicate_count = 0
while i < len(all_species):
    species = all_species[i]
    wavelength = all_wavelengths[i]
    depth = all_depths[i]
    
    if species in all_species[0:i] and wavelength in all_wavelengths[0:i]:
        duplicate_count += 1
        index = all_wavelengths.index(wavelength)
        duplicate_depth = all_depths[index]
        i -= 1

        #If depth of current line is larger than an existing line, delete the existing line
        if duplicate_depth < depth:
            duplicates_file.write(str(all_species[index]) + " " + str(all_wavelengths[index]) + " " + str(all_depths[index]) + "\n")
            del all_species[index]
            del all_wavelengths[index]
            del all_elements[index]
            del all_ionizations[index]
            del all_ionization_stages[index]
            del all_excitation_potentials[index]
            del all_log_gfs[index]
            del all_term_lowers[index]
            del all_term_uppers[index]
            del all_collisions[index]
            del all_gamma_radiatives[index]
            del all_gamma_starks[index]
            del all_gamma_van_der_waals[index]
            del all_widths[index]
            del all_depths[index]
            del all_lande_factors[index]

            del spec_ions[index]
            del wavelengths[index]
            del log_gfs[index]
            del E_lows[index]
            del J_lows[index]
            del E_ups[index]
            del J_ups[index]
            del lower_landes[index]
            del upper_landes[index]
            del mean_landes[index]
            del rads[index]
            del starks[index]
            del waalses[index]
            del depths[index]
            del line2s[index]
            del line3s[index]
            del line4s[index]            

        #If depth of current line is smaller than an existing line, delete the current line
        else:
            duplicates_file.write(str(all_species[i]) + " " + str(all_wavelengths[i]) + " " + str(all_depths[i]) + "\n")
            del all_species[i]
            del all_wavelengths[i]
            del all_elements[i]
            del all_ionizations[i]
            del all_ionization_stages[i]            
            del all_excitation_potentials[i]
            del all_log_gfs[i]
            del all_term_lowers[i]
            del all_term_uppers[i]
            del all_collisions[i]            
            del all_gamma_radiatives[i]
            del all_gamma_starks[i]
            del all_gamma_van_der_waals[i]            
            del all_widths[i]
            del all_depths[i]
            del all_lande_factors[i]

            del spec_ions[i]
            del wavelengths[i]
            del log_gfs[i]
            del E_lows[i]
            del J_lows[i]
            del E_ups[i]
            del J_ups[i]
            del lower_landes[i]
            del upper_landes[i]
            del mean_landes[i]
            del rads[i]
            del starks[i]
            del waalses[i]
            del depths[i]
            del line2s[i]
            del line3s[i]
            del line4s[i]

    i += 1

#Write all VALD parameters into a table
for i in range(0, len(all_species)):
    TABLE.write("#" + str(spec_ions[i]) + " " + str(all_ionization_stages[i]))
    TABLE.write("#") 
    TABLE.write(wavelengths[i])
    TABLE.write("#")
    TABLE.write(log_gfs[i])
    TABLE.write("#")
    TABLE.write(E_lows[i])
    TABLE.write("#")
    TABLE.write(J_lows[i])
    TABLE.write("#")
    TABLE.write(E_ups[i])
    TABLE.write("#")
    TABLE.write(J_ups[i])
    TABLE.write("#")
    TABLE.write(lower_landes[i])
    TABLE.write("#")
    TABLE.write(upper_landes[i])
    TABLE.write("#")
    TABLE.write(mean_landes[i])
    TABLE.write("#")
    TABLE.write(rads[i])
    TABLE.write("#")
    TABLE.write(starks[i])
    TABLE.write("#")
    TABLE.write(waalses[i])
    TABLE.write("#")
    TABLE.write(depths[i])
    TABLE.write("#")
    TABLE.write(line2s[i])
    TABLE.write("#")
    TABLE.write(line3s[i])
    TABLE.write("#")
    TABLE.write(line4s[i])
    TABLE.write("#\n")

#Write all parameters to text file
for i in range(0, len(all_species)):
    species = all_species[i]
    wavelength = all_wavelengths[i]
    element = all_elements[i]
    ionization = all_ionizations[i]
    ionization_num = all_ionization_stages[i]
    log_gf = all_log_gfs[i]
    excitation_potential = all_excitation_potentials[i]
    term_lower = all_term_lowers[i]
    term_upper = all_term_uppers[i]
    collisions = all_collisions[i]
    rad = all_gamma_radiatives[i]
    stark = all_gamma_starks[i]
    waals = all_gamma_van_der_waals[i]
    width = all_widths[i]
    depth = all_depths[i]
    lande_factor = all_lande_factors[i]
    LINES_file.write('[' + str(element) + str(ionization) + ' ' + str(wavelength) + ']\n')
    LINES_file.write('   Element=' + str(element) + '\n')
    LINES_file.write('   Ionization stage=' + str(ionization_num) + '\n')
    LINES_file.write('   Wavelength=' + str(wavelength) + '\n')
    LINES_file.write('   Excitation potential= ' + str(excitation_potential) + ' ev\n')
    LINES_file.write('   Log(gf)=' + str(log_gf) + '\n')
    LINES_file.write('   Term (lower)=' + str(term_lower) + '\n')
    LINES_file.write('   Term (upper)=' + str(term_upper) + '\n')
    LINES_file.write('   Collisions=1\n')
    LINES_file.write('   Gamma radiative=' + str(rad) + '\n')
    LINES_file.write('   Gamma stark=' + str(stark) + '\n')
    LINES_file.write('   Gamma van der Waals=' + str(waals) + '\n')
    LINES_file.write('   width=2\n')
    LINES_file.write('#Lande Factor = ' + str(lande_factor) + '\n')
    LINES_file.write('\n')

print('CONVERSION COMPLETE!')
print('TOTAL LINE COUNT: ' + str(total_line_count))
print('TOTAL ERROR COUNT: ' + str(total_error_count))
print('TOTAL DUPLICATE COUNT: ' + str(duplicate_count))

LINES_file.close()
error_file.close()
TABLE.close()
duplicates_file.close()
