import csv
import os

def transform_csv(input_filepath, output_filepath):
    # Open the original CSV file
    with open(input_filepath, 'r') as infile:
        reader = csv.reader(infile)
        
        # Extract the single row of values
        row = next(reader)
        
        # Split the row into individual values if it's a single string
        values = row[0].split(',') if len(row) == 1 else row
    
    # Open a new CSV file to write each value on a new line
    with open(output_filepath, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        
        # Write each value on a new row
        for value in values:
            writer.writerow([value])

# Define the directory containing the input CSV files
input_dir = '../../matRad-master/'
output_dir = './'  # Output files will be in the current directory

# List of filenames to process
filenames = ['coordsX.csv', 'coordsY.csv', 'coordsZ.csv']

# Process each file
for filename in filenames:
    input_filepath = os.path.join(input_dir, filename)
    output_filepath = os.path.join(output_dir, filename.replace('.csv', '_transformed.csv'))
    
    transform_csv(input_filepath, output_filepath)
