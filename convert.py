import sys

'''The following script changes the magnetic field values of a NICOLE model atmosphere file. 
Columns correspond to (according to NICOLE manual pg 30) log(tau_5000), 
temperature (in K), electron pressure (in dyn/cm^2), microturbulence 
(in cms), longitudinal magnetic field (gauss), line of sight velocity 
(cm/s) transverse field in x (gauss) and transverse field in y (gauss)'''

'''To run the following script, run "python convert.py <Input Model> <Output Model> x=<Bx> y=<By> z=<Bz>, where <Input Model> is the model you want to change, <Output Model> is the model you want to create, and Bx, By, and Bz are the magnetic fields (in G) in the x, y, and z directions, respectively. It is possible to change other columns as well, using tau, temp, ep, mt, or losv'''


columns = ['tau', 'temp', 'ep', 'mt', 'z', 'losv', 'x', 'y']
values = ['df']*8

#Get columns to change from terminal arguements and their corresponding values
for i in range(3, len(sys.argv)):
    arg = sys.argv[i]
    column = arg.split("=")[0]
    value = arg.split("=")[1]
    index = columns.index(column)
    values[index] = value

# Keep list of values to change
changes = []
for j in range(0, len(values)):
    if (values[j] != 'df'):
        changes.append(j)

# Open files given as arguements
f = open(str(sys.argv[1]), "r")
output = open(str(sys.argv[2]), "w+")
Lines = f.readlines()

#Iterate through input file, make necesary changes, and write to output file
for i in Lines:
    words = i.split()
    if(len(words) < 8):
        output.write(i)
    else:
        for i in range(0, len(values)):
            if values[i] != "df":
                words[i] = values[i]
        newline = " "
        newline = newline.join(words)
        output.write("  " + str(newline) + '\n')

output.close()


