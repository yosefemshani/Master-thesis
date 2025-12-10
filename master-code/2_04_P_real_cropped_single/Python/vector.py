import math

# Define pairs of vectors as tuples
pairs = [
    ((22.19, 122.87), (21.26, 122.89)),
    ((76.83, 125.79), (76.11, 125.87)),
    ((156.63, 134.13), (156.38, 133.97)),
    ((256.58, 150.62), (256.21, 149.24))
]

# Calculate magnitude and increment for each pair
results = []
for v1, v2 in pairs:
    # Magnitude of v1 and v2
    mag_v1 = math.sqrt(v1[0]**2 + v1[1]**2)
    mag_v2 = math.sqrt(v2[0]**2 + v2[1]**2)
    
    # Increment (Delta v)
    delta_v = mag_v1 - mag_v2
    
    # Store results
    results.append((mag_v1, mag_v2, delta_v))

print(results)