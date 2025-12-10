import numpy as np
import matplotlib.pyplot as plt

# Load the dose distribution along the beam path from the CSV file
dose_along_z = np.loadtxt('../../matRad-master/dose_along_beam_bixel_88.csv', delimiter=',')

# Plot the dose along the z-axis (depth) to visualize the Bragg peak
plt.figure(figsize=(8, 6))
plt.plot(dose_along_z, '-o', color='red')
plt.title('Bragg Peak (Dose Distribution Along Beam Path for Bixel 88)')
plt.xlabel('Depth (z-axis slice)')
plt.ylabel('Dose')
plt.grid(True)
plt.show()
