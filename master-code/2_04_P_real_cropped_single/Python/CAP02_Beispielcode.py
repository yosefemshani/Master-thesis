# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 15:01:46 2023

@author: Yosef
"""



import pandas as pd # to work with structured, tabulated data (.csv, .xlsx,..)
import numpy as np # to work fast with matrices, arrays and much more
from matplotlib import pyplot as plt # to visualize/plot your data

# read in the csv-file
outputfile_proton = '../TOPAS/TOPAS_Beam.csv'

# looking at the original csv-file, we know that the TOPAS simulation details
# are in first few lines (header) of the data sheet. We read in the first 8 
# lines using pandas to read the csv-file and setting nrows=7, since Python 
# starts counting with 0 instead of 1.
header = pd.read_csv(outputfile_proton, nrows = 7)

# Use the pandas package to read in the CSVfile as DataFrame. Skip all header
# lines, which all begin with #. Now, there is only data in the file, 
# therefore header=None.
df = pd.read_csv(outputfile_proton, comment='#', header=None)

#print(df)


proton_datamatrix = np.array(df)# convert dataframe df to array 

xrange = proton_datamatrix[:,0]
yrange = proton_datamatrix[:,1] #* 0.25
dose = proton_datamatrix[:,3]

# Filter the points where dose is non-zero
nonzero_dose_indices = np.where(dose > 0)
x_nonzero = xrange[nonzero_dose_indices]
y_nonzero = yrange[nonzero_dose_indices]

# Create a scatter plot to show where the dose is deposited
plt.figure(figsize=(10, 8))
plt.scatter(x_nonzero, y_nonzero, s=1, color='blue') # s=1 for small points
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Dose Deposition Locations')
plt.xlim(0, 340)
plt.ylim(0, 224)
plt.gca().invert_yaxis() # Invert y-axis to match typical image coordinates
plt.show()
plt.close()

#plt.figure(1) # creates a figure in which we plot 
#plt.plot(depth, dose, ls = '-', color = 'b') # plot depth on x, dose on y 
#plt.title('Tiefendosiskurve von Protonen mit 90 MeV.') # Title
#plt.xlabel('Tiefe / cm') # label for x-axis plt.ylabel('? / ?') # label for y axis plt.show()# show me the figure
#plt.ylabel('Dosis / Gy')
#plt.savefig("topasbeam.pdf")
#plt.close()

#dose_max = np.max(dose) # find dose maximum
#dose_norm = dose/dose_max 
#
#plt.figure(1) # creates a figure in which we plot 
#plt.plot(depth, dose_norm, ls = '-', color = 'b', label='TDK Proton') # plot depth on x, dose on y 
#plt.title('Tiefendosiskurve von 100 Protonen') # Title
#plt.xlabel('Tiefe / cm') # label for x-axis plt.ylabel('? / ?') # label for y axis plt.show()# show me the figure
#plt.ylabel('relative Dosis')
#plt.legend()
#plt.savefig("proton100.pdf")
#plt.close()