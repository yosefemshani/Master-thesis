# pixel2mm.py

import csv

def convert_pixel_to_mm(input_file, output_file, scale_factor):
    with open(input_file, mode='r') as infile:
        reader = csv.reader(infile)

        # Open the output file to write the modified data
        with open(output_file, mode='w', newline='') as outfile:
            writer = csv.writer(outfile)

            # Iterate through each row in the input file (start from the first entry)
            for row in reader:
                x = float(row[0]) * scale_factor  # Scale the x value
                y = float(row[1]) * scale_factor  # Scale the y value
                energy = float(row[2])  # Leave energy value unchanged

                # Write the updated values to the output file
                writer.writerow([x, y, energy])

if __name__ == "__main__":
    input_filename = "../../matRad-master/maxRows_3D_with_energy.csv"
    output_filename = "../../matRad-master/maxRows_mm_with_energy.csv"
    scale_factor = 1.09375

    convert_pixel_to_mm(input_filename, output_filename, scale_factor)
