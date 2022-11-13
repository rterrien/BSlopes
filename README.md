# BSlopes

## Features:

This tool allows users to generate magnetic sensitivities for a linelist that they input. In addition, it also provides some helpful intermediary tools to convert, for example, a VALD output to a NICOLE LINES file or MARCS model atmosphere to a NICOLE model atmosphere. 

## Dependencies:

The tool makes use of the following:

- NICOLE spectral synthesis code: https://github.com/hsocasnavarro/NICOLE
- Pandas, PyAstronomy, Scipy.stats, NumPy, Matplotlib

## Overview:

The tool contains the following Python scripts:

- VALD_to_NICOLE.py: converts VALD outputs (using the extract stellar option, long format) to NICOLE LINES files and also outputs textfiles of duplicate lines, lines with errors, and VALD parameters in a table format called “VALD_TABLE.txt”
- create_model_atmospheres.py: creates a list of NICOLE model atmospheres based on a list of temperatures (and corresponding MARCS model atmospheres) and magnetic field inputs
- MARCS_to_NICOLE.py: converts a MARCS model atmosphere into a NICOLE model atmosphere
- convert.py: changes a NICOLE model atmosphere by changing the magnetic field values
- run_VALD_lines.py: takes a NICOLE LINES file and NICOLE model atmospheres as input and runs NICOLE on all lines
- line_fitting.py: reads a list of model outputs from NICOLE, determines changes in EW, FWHM, and Depth by varying the magnetic field and temperature, and writes these slopes to a text file called “SLOPES.txt”
- create_table.py: takes as input VALD parameters in a text file called “VALD_TABLE.txt” and magnetic sensitivities in a text file called “SLOPES.txt” and combines the tables to create “FINAL_TABLE.csv”

## Steps:

1. Get all the dependencies
2. In your NICOLE folder, navigate to the test directory and copy the “syn1” directory. 
3. Clone this repository in the copy of the syn1 directory
4. In the INPUTS directory, add the MARCS model atmospheres 
5. In the INPUTS directory, change the following three input text files:
    - magnetic_field_input.txt: a list of B (in G) for the simulation
    - temperature_input.txt: a list of T (in K) and the corresponding MARCS atmosphere file names that were added in step 4
    - VALD_input.txt: the output of a VALD query using the long format
6. Run ```python simulate.py``` in the copy of the syn1 directory

## Outputs:

The tool creates a directory called OUTPUT, followed by the date and time when simulate.py was run. In this OUTPUT directory are the following:

- FINAL_TABLE.csv, containing a list of lines, their VALD parameters, their EW, FWHM, ad depth values at various T and B, and their T and B sensitivities
- A Figures directory, the varying T and B line profiles for all lines, the line fits of EW, FWHM, and Depth vs T and B for all lines, and scatter plots of all EW, FWHM, and Depth sensitivities vs the Landé Factor. 
- An INPUTS directory containing the INPUTS used
- Additional files that were used in the simulation, such as model_atmospheres, model_files, NICOLE LINES file, SLOPES.txt, DUPLICATES.txt, ERRORS.txt, and VALD_TABLE.txt


