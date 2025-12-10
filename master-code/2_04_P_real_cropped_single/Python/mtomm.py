import csv

# Function to convert centimeters to millimeters
def centimeters_to_millimeters(value):
    return value * 10  # 1 centimeter = 10 millimeters

# Input and output file paths
input_file = '../../matRad-master/trajectory.csv'
output_file = 'trajectory_mm.csv'

# Read input CSV file and write output CSV file with values converted to millimeters
with open(input_file, 'r') as csv_in, open(output_file, 'w', newline='') as csv_out:
    reader = csv.reader(csv_in)
    writer = csv.writer(csv_out)
    
    for row in reader:
        # Convert values from centimeters to millimeters
        x_mm = centimeters_to_millimeters(float(row[0]))
        y_mm = centimeters_to_millimeters(float(row[1]))
        
        # Write the converted values to the output CSV file
        writer.writerow([x_mm, y_mm])

print("Conversion completed. Output written to", output_file)
