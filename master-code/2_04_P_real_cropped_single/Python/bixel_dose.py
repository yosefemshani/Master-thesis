import numpy as np
import matplotlib.pyplot as plt

# Load the 2D dose distribution data for slice z=2 from the CSV file
dose_slice_z2 = np.loadtxt('../../matRad-master/physicalDoseData_slice_z_2_bixel_88.csv', delimiter=',')

# Plot the 2D dose distribution for slice z=2
plt.figure(figsize=(8, 6))
plt.imshow(dose_slice_z2, cmap='hot', origin='lower')
plt.colorbar(label='Dose')
plt.title('2D Dose Distribution (Slice z=2) for Bixel 88')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
