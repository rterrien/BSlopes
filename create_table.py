import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from PyAstronomy import pyasl
from scipy.stats import gaussian_kde
import subprocess

'''The following script conmbines a textfile of slopes and a texfile of VALD outputs. It also plots the slopes of all lines vs the Lande factor.'''

'''To run the following script, run "python VALD_to_NICOLE.py <Slopes> <VALD outputs>", where <Slopes> is a textfile listing the slopes of every line and <VALD outputs> lists the VALD output parameters used for every line.'''

#Read the textfile containing the slopes
df=pd.read_csv(sys.argv[1], sep=' ', error_bad_lines=False)
df = df.loc[:,~df.columns.str.match("Unnamed")]

#Read the textfile containing VALD outputs
data = pd.read_csv(sys.argv[2], delimiter="#", header=None)

#Clean the columns
data2 = data.iloc[::4, :]
data2.columns = ['Element', 'Ionization', 'Wavelength (Air, Å)', 'log gf', 'E low (eV)', 'J low', 'E up (eV)', 'J up', 'Lande lower', 'Lande upper', 'Lande mean', 'Radiative', 'Stark', 'Waals', 'Depth', 'Lower term']
data2[['Element', 'Ionization']] = data2['Ionization'].str.split(' ', 1, expand=True)
sources = data2['Lower term']
sources_stripped = []
for i in sources:
    sources_stripped.append(' '.join(i.split()))
data2['Lower term'] = sources_stripped
data3 = data.iloc[1::4, 1]
data4 = data.iloc[2::4, 1]
data3_stripped = []
for i in data3:
    data3_stripped.append(' '.join(i.split()))
data4_stripped = []
for i in data4:
    data4_stripped.append(' '.join(i.split()))
data2['Upper term'] = data3_stripped
data2['Source'] = data4_stripped

#Set the wavelenth as the index for the VALD dataframe, and round values
data2.set_index('Wavelength (Air, Å)', inplace=True, drop = False)
VALD_index = data2.index
rounded = [np.round(x, 10) for x in VALD_index]
data2['Wavelength (Air, Å)'] = rounded
data2.set_index('Wavelength (Air, Å)', inplace=True, drop = False)

#Set the wavelenth as the index for the slopes dataframe, and round values
Slopes_index = df['Wavelength']
rounded = [np.round(x, 10) for x in Slopes_index]
df['Wavelength'] = rounded
df.set_index('Wavelength', inplace=True, drop = False)

#Merge the VALD and slope datafiles, using the index of both dataframes
TABLE = pd.merge(data2, df, left_index=True, right_index=True, how = 'inner')
TABLE.drop(columns = ['Species', 'Wavelength', 'VALD_Lande_Factor'], inplace=True)

#Save the final table as a .CSV
TABLE.to_csv('OUTPUTS/FINAL_TABLE.csv', index=False)

#Define function that colors points in a scatter plot by the density at that location
def get_density(xx, yy):
    x = []
    y = []
    for i in range(0, len(xx)):
        x.append(xx[i]/np.mean(xx))
        y.append(yy[i])
    x = np.array(x)
    x = x.flatten()
    y = np.array(y)
    y = y.flatten()
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    return x, y, z

#Get the Lande factors
Lande_mean_index = TABLE.columns.get_loc('Lande mean')
Lande_estimated_index = TABLE.columns.get_loc('Estimated_Lande_Factor')

#Get a list of the slope columns
slope_lists = []
for i in TABLE.columns:
    if 'vs' in i:
        slope_lists.append(i)
lande_factors = []

#If the VALD Lande factor is not physical (within the arbitrary range -1 to 4), use the estimated Lande factor
for j in range(0, len(TABLE)):
    if TABLE.iloc[j, Lande_mean_index] > 4 or TABLE.iloc[j, Lande_mean_index] < -1:
        lande_factors.append(TABLE.iloc[j, Lande_estimated_index])
    else:
        lande_factors.append(TABLE.iloc[j, Lande_mean_index])

#Create a directory that will contain all the slopes vs Lande Factor figures
subprocess.call(['mkdir', 'OUTPUTS/Figures/SlopesvsLande'])

#For each slope, plot slopes of all lines vs the Lande factor
figure = plt.figure(figsize = (15, 10))
for i in slope_lists:
    plt.cla()
    plt.clf()
    if len(TABLE[i].tolist()) > 1:
        x, y, z = get_density(lande_factors, TABLE[i].tolist())
        plt.scatter(x, y, c = z, s = 10)
    else:
        plt.scatter(lande_factors, TABLE[i].tolist())
    plt.xlabel("Lande Factor", size = 14)
    plt.ylabel(i, size = 14)
    plt.xlim(-1, 3)
    fig = plt.gcf()
    fig.set_size_inches(15, 10)
    plt.tight_layout()
    plt.savefig("OUTPUTS/Figures/SlopesvsLande/" + str(i) + ".png")

