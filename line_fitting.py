from __future__ import division 
import model_prof_tools
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm     as cm
from PyAstronomy import pyasl
#import hpfspec2
import sys
from scipy import stats
from matplotlib.ticker import FormatStrFormatter 

'''The following script iterates through all .pro files and finds the EW, FWHM, and Depth of each line profile. It also performs a linear fit to each parameter against the temperature and magnetic field and writes the results to a text file. Plots of varying temperature and magnetic field are also created, as well as the linear fits to the structural parameters of the line agagainst temperature and magnetic field.'''

'''To run the following script, run "python line_fitting.py <Temperature input file> <Magnetic Field input file> <LINES>" where <Temperature input file> is a text file, with each line containing a temperature value and a MARCS model atmosphere (separated with a space), <Magnetic Field input file> is a text file, with each line containing a magnetic field value, and <LINES> is a textfile containing a list of LINES.'''

T_input = open(sys.argv[1], "r")
T_lines = T_input.readlines()
B_input = open(sys.argv[2], "r")
B_lines = B_input.readlines()
temps = []
mfs = []
MARCS_model_atmosphere_list = []

#Read tempature values and MARCS model atmospheres from the temperature input file
for i in T_lines:
    words = i.split()
    temps.append(int(words[0]))
    MARCS_model_atmosphere_list.append(words[1])

#Read magnetic field values from the magnetic field input file
for j in B_lines:
    words = j.split()
    mfs.append(int(words[0]))

#The parameters to parse through the wavelengths
lambda_step = '2'
num_lambda = '2000'

#Saves a figure depicting the change of the given line with varying magnetic field
def plot_all_mf(element, central_wavelength):
    plt.cla()
    plt.clf()
    colormap = cm.plasma(np.linspace(0, 1, len(mfs)))
    for i in range(0, len(mfs)):
        ws, ints = get_wavelength_intensity_mf(i)
        plt.plot(ws, ints, label = str(mfs[i]) + " G", color = colormap[i], linewidth = 4)
    sm = plt.cm.ScalarMappable(cmap="plasma", norm=plt.Normalize(vmin=min(mfs), vmax=max(mfs)))
    cbar = plt.colorbar(sm)
    cbar.set_label('Magnetic Field Strength (G)', size = 20)
    cbar.ax.tick_params(labelsize=14)
    #Fit delta lambda such that at 3700 angstroms, the offset is 0.2 angstroms and at 15000 ansgtroms, the offset is 2 angstroms
    ws_offset = 1.59e-4 * float(central_wavelength) - 0.3883
    plt.xlim(float(central_wavelength) - ws_offset, float(central_wavelength) + ws_offset)
    plt.title(str(element) + " " + str(round(float(central_wavelength), 4)), size = 20)
    plt.xlabel(r"Wavelength $(\AA)$", size = 20)
    plt.ylabel("I / $I_c$", size = 20)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    fig = plt.gcf()
    fig.set_size_inches(15, 10)
    plt.tight_layout()
    plt.savefig("OUTPUTS/Figures/varying_B/" + str(element) + "_" + str(central_wavelength) + ".png")

#Saves a figure depicting the change of the given line with temperature
def plot_all_temp(element, central_wavelength):
    plt.cla()
    plt.clf()
    colormap = cm.plasma(np.linspace(0, 1, len(temps)))
    for i in range(0, len(temps)):
        ws, ints = get_wavelength_intensity_temp(i)
        plt.plot(ws, ints, label = "T = " + str(int(i)) + " K", color = colormap[i], linewidth = 4)
    sm = plt.cm.ScalarMappable(cmap="plasma", norm=plt.Normalize(vmin=min(temps), vmax=max(temps)))
    cbar = plt.colorbar(sm)
    cbar.set_label('Temperature (K)', size = 20)
    cbar.ax.tick_params(labelsize=14) 

    #Fit delta lambda such that at 3700 angstroms, the offset is 0.2 angstroms and at 15000 ansgtroms, the offset is 2 angstroms
    ws_offset = 1.59e-4 * float(central_wavelength) - 0.3883

    plt.xlim(float(central_wavelength) - ws_offset, float(central_wavelength) + ws_offset)
    plt.title(str(element) + " " + str(round(float(central_wavelength), 4)), size = 20)
    plt.xlabel(r"Wavelength $(\AA)$", size = 20)
    plt.ylabel("I / $I_c$", size = 20)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    fig = plt.gcf()
    fig.set_size_inches(15, 10)
    plt.tight_layout()
    plt.savefig("OUTPUTS/Figures/varying_T/" + str(element) + "_" + str(central_wavelength) + ".png")

#Gets the wavelength and intensity values of a given line (with varying magnetic field)
def get_wavelength_intensity_mf(num):
    wavelength, intensity = loader(mf_file_names[num])
    return wavelength, intensity

#Gets the wavelength and intensity values of a given line (with varying magnetic field)
def get_wavelength_intensity_temp(num):
    wavelength, intensity = loader(temp_file_names[num])
    return wavelength, intensity

#Helper function to load a file name and parse the different parameters in the name
def loader(file_name):
    words = file_name.replace('OUTPUTS/model_files/', '').split("_")
    wstart = float(words[1]) - 2
    wstep = int(words[4])/1000
    nsteps = int(words[5])
    ws = np.arange(wstart,wstart+wstep*nsteps,wstep)
    file_name_2 = words[0] + '_' + words[1].rstrip("0") + '_' + words[2] + '_' + words[3] + '_' + words[4] + '_' + words[5]
    inp = model_prof_tools.read_prof('OUTPUTS/model_files/' + file_name_2 + '.pro','nicole',0,0,nsteps,0,0)
    inp2 = np.reshape(inp,(nsteps,4))
    intensity = inp2[:,0]/np.nanmax(inp2[:,0])
    return ws, intensity

#Gets the line profile for a particular temperature
def get_line_for_temp(temp):
    if temp < min(temps):
        print("temperature too low")
    elif temp > max(temps):
        print("temperature too high")
    else:
        line_low = 0
        line_high = 0
        found = False
        for i in range(0, len(temps)):
            if temp == temps[i] and found == False:
                found = True
                ws, ints = get_wavelength_intensity_temp(i)
                return ws, ints
            elif temp <= temps[i] and found == False:
                found = True
                line_low = i - 1
                line_high = i
                ws_low, ints_low = get_wavelength_intensity_temp(line_low)
                ws_high, ints_high = get_wavelength_intensity_temp(line_high)
                low_diff = temp - temps[line_low]
                high_diff = temps[line_high] - temp
                low_r = low_diff / (low_diff + high_diff)
                high_r = high_diff / (low_diff + high_diff)
                index = 150
                ws = (ws_low * (1 - low_r) + ws_high * (1 - high_r))
                ints = (ints_low * (1 - low_r) + ints_high * (1 - high_r))
                return ws, ints 
            
#Gets the line profile for a particular magnetic field
def get_line_for_mf(mf):
    index = mfs.index(mf)
    ws, ints = get_wavelength_intensity_mf(index)
    return ws, ints

#Gets the equivalent width of a line, given wavelength, intensities, and left and right boundaries of the line
def get_equivalent_width(wl,fl,limit_left,limit_right):
    assert limit_left > np.nanmin(wl)
    assert limit_right < np.nanmax(wl)
    bin_size = np.diff(wl)
    bin_size = np.concatenate(([bin_size[0]],bin_size))
    bin_left = wl - bin_size/2.
    bin_right = wl + bin_size/2.
    
    condition_finite = np.isfinite(fl)
    condition_all_in = (bin_left >= limit_left) & (bin_right <= limit_right)
    
    use = np.nonzero(condition_finite & condition_all_in)[0]
    wluse = wl[use]
    fluse = fl[use]
    
    bins = np.diff(wluse)
    bins = np.concatenate(([bins[0]],bins))
    
    sub = (1. - fluse) * bins
    
    leftmost_index = use[0]
    left_extra_bin = bin_right[leftmost_index-1] - limit_left
    left_extra_val = (1. - fl[leftmost_index-1]) * left_extra_bin
    
    rightmost_index = use[-1]
    right_extra_bin = limit_right - bin_left[rightmost_index+1]
    right_extra_val = (1. - fl[rightmost_index+1]) * right_extra_bin
    
    return(np.sum(sub) + left_extra_val + right_extra_val)

#Gets the full-width at half max of a line, given the wavelength and intensities of the line
def get_fwhm(ws, fluxes):
    mins = local_min(fluxes)
    mins2 = []
    threshold = 0.99
    for i in mins:
        if i < threshold:
            mins2.append(i)
    lambdas = []
    for i in mins2:
        index = fluxes.tolist().index(i)
        lambdas.append(ws[index])
    try:
        low_half_flux = 1 - (1 - mins2[0])/2
        high_half_flux = 1 - (1 - mins2[-1])/2
        half_flux = np.array([low_half_flux]*len(ws))  
        xs, ys = interpolated_intercepts(ws, fluxes, half_flux)
        lower_x, lower_y = xs[0], ys[0]
        higher_x, higher_y = xs[-1], ys[-1]
        fwhm = higher_x - lower_x
        return(fwhm[0])
    except:
        return 0

#Gets all local minima of a line
def local_min(ys):
    return [y for i, y in enumerate(ys)
            if ((i == 0) or (ys[i - 1] >= y))
            and ((i == len(ys) - 1) or (y < ys[i+1]))]

#Helper function for the get_fwhm() function to interpolate two curves and fine the approximate intercept, see https://stackoverflow.com/questions/42464334/find-the-intersection-of-two-curves-given-by-x-y-data-with-high-precision-in
def interpolated_intercepts(x, y1, y2):
    """Find the intercepts of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """    

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idxs = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)

    xcs = []
    ycs = []
    

    for idx in idxs:
        xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
        xcs.append(xc)
        ycs.append(yc)
    return np.array(xcs), np.array(ycs)

#Gets the depth of a line based on the wavelength and intensities of the line
def get_line_depth(ws, fluxes):
    minimum = min(fluxes)
    return 1 - minimum

#Read the LINES file
VALD_LINES_file = open(sys.argv[3], "r")
VALD_LINES = VALD_LINES_file.readlines()
species = []
wavelengths = []
lande_factors = []
lande_factors_estimated = []
L_List = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N']
            
#Estimate the Lande Factor using Kochukhov 2021       
count = 0
for i in range(0, len(VALD_LINES)):
    line = VALD_LINES[i]
    words = line.split()
    if len(words) > 0 and words[0][0] == '[' and words[-1][-1] == ']':
        species.append(words[0][1::])
        wavelengths.append(words[-1][0:-1])
    elif len(words) > 1 and words[0] == '#Lande':
        lande_factors.append(float(words[3]))
    elif len(words) > 1 and words[0] == 'Term':
        if words[1].split('=')[0] == '(lower)':
            value = words[1].split('=')[1]
            L_lower = float(L_List.index(value[1]))
            S_lower = (float(value[0]) - 1)/2
            J_lower = float(value[2::])
        elif words[1].split('=')[0] == '(upper)':
            value = words[1].split('=')[1]
            L_upper = float(L_List.index(value[1]))
            S_upper = (float(value[0]) - 1)/2
            J_upper = float(value[2::])
            if J_lower == 0 and J_upper != 0:
                g_lower = 0
                g_upper = 3/2 + (S_upper*(S_upper + 1) - L_upper*(L_upper + 1))/(2*J_upper*(J_upper + 1))
            elif J_upper == 0 and J_lower != 0:
                g_upper = 0
                g_lower = 3/2 + (S_lower*(S_lower + 1) - L_lower*(L_lower + 1))/(2*J_lower*(J_lower + 1))
            elif J_lower == 0 and J_upper == 0:
                g_lower = 0
                g_upper = 0
            else:
                g_lower = 3/2 + (S_lower*(S_lower + 1) - L_lower*(L_lower + 1))/(2*J_lower*(J_lower + 1))
                g_upper = 3/2 + (S_upper*(S_upper + 1) - L_upper*(L_upper + 1))/(2*J_upper*(J_upper + 1))

            geff = 1/2 * (g_lower + g_upper) + 1/4 * (g_lower - g_upper)*(J_lower*(J_lower + 1) - J_upper*(J_upper + 1))
            lande_factors_estimated.append(float(geff))
            count += 1

#Open file to write the slopes
subprocess.call(['touch', 'OUTPUTS/SLOPES.txt'])
file_object = open("OUTPUTS/SLOPES.txt", "w")
columns = ["Species", "Wavelength", "VALD_Lande_Factor", "Estimated_Lande_Factor"]

#Print the column names
for i in range(0, len(temps)):
    columns.append('EWT_' + str(int(temps[i])))
    columns.append('FWHMT_' + str(int(temps[i])))
    columns.append('DepthT_' + str(int(temps[i])))

for i in range(0, len(mfs)):
    columns.append('EWB_' + str(int(mfs[i])))
    columns.append('FWHMB_' + str(int(mfs[i])))
    columns.append('DepthB_' + str(int(mfs[i])))

for i in range(0, len(temps) - 1):
    columns.append('EWvsT_' + str(int(temps[i])) + "_" + str(int(temps[i+1])))
    
for i in range(0, len(mfs) - 1):
    columns.append('EWvsB_' + str(int(mfs[i])) + "_" + str(int(mfs[i+1])))
    
for i in range(0, len(temps) - 1):
    columns.append('FWHMvsT_' + str(int(temps[i])) + "_" + str(int(temps[i+1])))
    
for i in range(0, len(mfs) - 1):
    columns.append('FWHMvsB_' + str(int(mfs[i])) + "_" + str(int(mfs[i+1])))
    
for i in range(0, len(temps) - 1):
    columns.append('DepthvsT_' + str(int(temps[i])) + "_" + str(int(temps[i+1])))
    
for i in range(0, len(mfs) - 1):
    columns.append('DepthvsB_' + str(int(mfs[i])) + "_" + str(int(mfs[i+1])))
    
for i in range(0, len(temps) - 1):
    columns.append('EWfvsT_' + str(int(temps[i])) + "_" + str(int(temps[i+1])))
    
for i in range(0, len(mfs) - 1):
    columns.append('EWfvsB_' + str(int(mfs[i])) + "_" + str(int(mfs[i+1])))
    
for i in range(0, len(temps) - 1):
    columns.append('FWHMfvsT_' + str(int(temps[i])) + "_" + str(int(temps[i+1])))
    
for i in range(0, len(mfs) - 1):
    columns.append('FWHMfvsB_' + str(int(mfs[i])) + "_" + str(int(mfs[i+1])))
    
for i in range(0, len(temps) - 1):
    columns.append('DepthfvsT_' + str(int(temps[i])) + "_" + str(int(temps[i+1])))
    
for i in range(0, len(mfs) - 1):
    columns.append('DepthfvsB_' + str(int(mfs[i])) + "_" + str(int(mfs[i+1])))

for i in columns:
    file_object.write(str(i) + " ")
file_object.write('\n')

#Iterate through all the lines and get their EW, FWHM, and Depth
FULL_RANGE = np.arange(0, len(species)).tolist()
bad_lines2 = []
figure = plt.figure(figsize = (15, 10))

for k in FULL_RANGE:

    element = species[k]
    central_wavelength = wavelengths[k]
    try:
        central_wavelength_vac = pyasl.airtovac2(central_wavelength)
    except:
        #central_wavelength_vac = pyasl.airtovac2(central_wavelength, precision = 1e-10)
        continue
        
    lande_factor = lande_factors[k]
    lande_factor_estimated = lande_factors_estimated[k]
    temp_file_names = []
    mf_file_names = []
    
    #List all the model files that vary with temperature
    for i in temps:
        file_name = 'OUTPUTS/model_files/' + str(element) + "_" + str(central_wavelength) + "_T" + str(i) + "_" + "MF0" + "_" + lambda_step + "_" + num_lambda
        temp_file_names.append(file_name)
    
    #List all the model files that vary with magnetic field
    for i in mfs:
        file_name = 'OUTPUTS/model_files/' + str(element) + "_" + str(central_wavelength) + "_T" + str(temps[0]) + "_MF" + str(i) +  "_" + lambda_step + "_" + num_lambda
        mf_file_names.append(file_name)
    
    try:
        passes_depth_check = True
        file_object.write(str(element) + ' ' + str(central_wavelength) + ' ')
        file_object.write(str(lande_factor) + ' ')
        file_object.write(str(round(lande_factor_estimated, 3)) + ' ')
        
        #Plot line profiles with varying temperature and magnetic field
        plot_all_temp(element, central_wavelength)
        plot_all_mf(element, central_wavelength)

        #Arbitrary wavelength offset for left and right boundaries for the EW calculation
        ws_offset = 1.95

        ew_temp = []
        fwhm_temp = []
        depth_temp = []

        #Get the EW, FWHM, and Depth for the current line for varying temperature
        for i in temps:
            ws, ints = get_line_for_temp(i)
            min_int = min(ints)
            index = ints.tolist().index(min_int)
            min_w = ws[index]
            ew = get_equivalent_width(ws, ints, float(central_wavelength) - ws_offset, float(central_wavelength) + ws_offset)
            ew_temp.append(ew)
            fwhm = get_fwhm(ws, ints)
            if fwhm == 0:
                passes_depth_check = False
            fwhm_temp.append(fwhm)
            line_depth = get_line_depth(ws, ints)
            depth_temp.append(line_depth)
            file_object.write(str(ew) + " ")
            file_object.write(str(fwhm) + " ")
            file_object.write(str(line_depth) + " ")
        
        ew_mf = []
        fwhm_mf = []
        depth_mf = []

        #Get the EW, FWHM, and Depth for the current line for varying magnetic field
        for i in mfs:
            ws, ints = get_line_for_mf(i)
            ew = get_equivalent_width(ws, ints, float(central_wavelength) - ws_offset, float(central_wavelength) + ws_offset)
            ew_mf.append(ew)
            fwhm = get_fwhm(ws, ints)
            if fwhm == 0:
                passes_depth_check = False
            fwhm_mf.append(fwhm)
            line_depth = get_line_depth(ws, ints)
            depth_mf.append(line_depth)
            file_object.write(str(ew) + " ")
            file_object.write(str(fwhm) + " ")
            file_object.write(str(line_depth) + " ")            
        
        if passes_depth_check:
           
            #Define fractional changes for the EW, FWHM, and Depth. Take the line profile simulated at the lowest temperature to be the control line  
            control_ew_temp = ew_temp[0]
            ew_temp_f = ew_temp - control_ew_temp 
            ew_temp_f = ew_temp_f / control_ew_temp
            control_fwhm_temp = fwhm_temp[0]
            fwhm_temp_f = fwhm_temp - control_fwhm_temp 
            fwhm_temp_f = fwhm_temp_f / control_fwhm_temp
            control_depth_temp = depth_temp[0]
            depth_temp_f = depth_temp - control_depth_temp 
            depth_temp_f = depth_temp_f / control_depth_temp
            ew_mf_f = ew_mf - control_ew_temp 
            ew_mf_f = ew_mf_f / control_ew_temp
            fwhm_mf_f = fwhm_mf - control_fwhm_temp 
            fwhm_mf_f = fwhm_mf_f / control_fwhm_temp
            depth_mf_f = depth_mf - control_depth_temp 
            depth_mf_f = depth_mf_f / control_depth_temp


            #Fit and plot the EW vs temperature
            plt.cla()
            plt.clf()            
            for i in range(0, len(temps) - 1):
                slope, intercept, r, p, se = stats.linregress(temps[i:i+2], ew_temp[i:i+2])
                file_object.write(str(slope) + ' ') 
            plt.scatter(temps, ew_temp)
            plt.plot(temps, ew_temp)
            plt.xlabel(r"Temperature (K)", size = 14)
            plt.ylabel("EW ($\AA$)", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/EWvsT/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()
            
            #Fit and plot the EW vs magnetic field
            plt.cla()
            plt.clf()
            for i in range(0, len(mfs) - 1):
                slope, intercept, r, p, se = stats.linregress(mfs[i:i+2], ew_mf[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(mfs, ew_mf)
            plt.plot(mfs, ew_mf)
            plt.xlabel(r"B (Gauss)", size = 14)
            plt.ylabel("EW ($\AA$)", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/EWvsB/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()

            #Fit and plot the FWHM vs temperature
            plt.cla()
            plt.clf()            
            for i in range(0, len(temps) - 1):
                slope, intercept, r, p, se = stats.linregress(temps[i:i+2], fwhm_temp[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(temps, fwhm_temp)
            plt.plot(temps, fwhm_temp)
            plt.xlabel(r"Temperature (K)", size = 14)
            plt.ylabel("FWHM ($\AA$)", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/FWHMvsT/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()
 
            #Fit and plot the FWHM vs magnetic field
            plt.cla()
            plt.clf()
            for i in range(0, len(mfs) - 1):
                slope, intercept, r, p, se = stats.linregress(mfs[i:i+2], fwhm_mf[i:i+2])
                #abline(slope, intercept)
                file_object.write(str(slope) + ' ')    
            plt.scatter(mfs, fwhm_mf)
            plt.plot(mfs, fwhm_mf)
            plt.xlabel(r"B (Gauss)", size = 14)
            plt.ylabel("FWHM ($\AA$)", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/FWHMvsB/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()
 
            #Fit and plot the Depth vs temperature
            plt.cla()
            plt.clf() 
            for i in range(0, len(temps) - 1):
                slope, intercept, r, p, se = stats.linregress(temps[i:i+2], depth_temp[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(temps, depth_temp)
            plt.plot(temps, depth_temp)
            plt.xlabel(r"Temperature (K)", size = 14)
            plt.ylabel("Depth ($\AA$)", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/DepthvsT/" + str(element) + "_" + str(central_wavelength) + ".png")             
            plt.close()

            #Fit and plot the Depth vs magnetic field
            plt.cla()
            plt.clf()
            for i in range(0, len(mfs) - 1):
                slope, intercept, r, p, se = stats.linregress(mfs[i:i+2], depth_mf[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(mfs, depth_mf)
            plt.plot(mfs, depth_mf)
            plt.xlabel(r"B (Gauss)", size = 14)
            plt.ylabel("Depth ($\AA$)", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/DepthvsB/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()

            #Fit and plot the EW vs temperature (fractional change)
            plt.cla()
            plt.clf()
            for i in range(0, len(temps) - 1):
                slope, intercept, r, p, se = stats.linregress(temps[i:i+2], ew_temp_f[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(temps, ew_temp_f)
            plt.plot(temps, ew_temp_f)
            plt.xlabel(r"Temperature (K)", size = 14)
            plt.ylabel("$\Delta$EW / EW$_c$", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/EWfvsT/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()

            #Fit and plot the EW vs magnetic field (fractional change)
            plt.cla()
            plt.clf()
            for i in range(0, len(mfs) - 1):
                slope, intercept, r, p, se = stats.linregress(mfs[i:i+2], ew_mf_f[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(mfs, ew_mf_f)
            plt.plot(mfs, ew_mf_f)
            plt.xlabel(r"B (Gauss)", size = 14)
            plt.ylabel("$\Delta$EW / EW$_c$", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/EWfvsB/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()
 
            #Fit and plot the FWHM vs temperature (fractional change)
            plt.cla()
            plt.clf()
            for i in range(0, len(temps) - 1):
                slope, intercept, r, p, se = stats.linregress(temps[i:i+2], fwhm_temp_f[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(temps, fwhm_temp_f)
            plt.plot(temps, fwhm_temp_f)
            plt.xlabel(r"Temperature (K)", size = 14)
            plt.ylabel("$\Delta$FWHM / FWHM$_c$", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/FWHMfvsT/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()

            #Fit and plot the FWHM vs magnetic field (fractional change)
            plt.cla()
            plt.clf()
            for i in range(0, len(mfs) - 1):
                slope, intercept, r, p, se = stats.linregress(mfs[i:i+2], fwhm_mf_f[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(mfs, fwhm_mf_f)
            plt.plot(mfs, fwhm_mf_f)
            plt.xlabel(r"B (Gauss)", size = 14)
            plt.ylabel("$\Delta$FWHM / FWHM$_c$", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/FWHMfvsB/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()

            #Fit and plot the Depth vs temperature (fractional change)
            plt.cla()
            plt.clf()
            for i in range(0, len(temps) - 1):
                slope, intercept, r, p, se = stats.linregress(temps[i:i+2], depth_temp_f[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(temps, depth_temp_f)
            plt.plot(temps, depth_temp_f)
            plt.xlabel(r"Temperature (K)", size = 14)
            plt.ylabel("$\Delta$Depth / Depth$_c$", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/DepthfvsT/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()

            #Fit and plot the Depth vs magnetic field (fractional change)
            plt.cla()
            plt.clf()
            for i in range(0, len(mfs) - 1):
                slope, intercept, r, p, se = stats.linregress(mfs[i:i+2], depth_mf_f[i:i+2])
                file_object.write(str(slope) + ' ')    
            plt.scatter(mfs, depth_mf_f)
            plt.plot(mfs, depth_mf_f)
            plt.xlabel(r"B (Gauss)", size = 14)
            plt.ylabel("$\Delta$Depth / Depth$_c$", size = 14)
            plt.title(str(element) + " " + str(central_wavelength))
            fig = plt.gcf()
            fig.set_size_inches(15, 10)
            plt.tight_layout()
            plt.savefig("OUTPUTS/Figures/DepthfvsB/" + str(element) + "_" + str(central_wavelength) + ".png")
            plt.close()
            file_object.write('\n')

    except:
        print("ERROR!")

file_object.close()
