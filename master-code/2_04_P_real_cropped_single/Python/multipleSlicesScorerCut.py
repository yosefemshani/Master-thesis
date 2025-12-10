# File path of the original and new CSV files
input_file = '../TOPAS/TOPAS_Beam.csv'
output_file = 'TOPAS_Beam_Slice1.csv'

# Number of lines to read (header + one slice of data)
lines_to_read = 76748

def extract_first_slice(input_file, output_file, lines_to_read):
    try:
        # Open the input file for reading
        with open(input_file, 'r') as infile:
            # Open the output file for writing
            with open(output_file, 'w') as outfile:
                # Read and write the specified number of lines
                for i in range(lines_to_read):
                    line = infile.readline()
                    if line:
                        outfile.write(line)
                    else:
                        break
        print(f"First slice extracted to {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

# Run the function to extract the first slice
extract_first_slice(input_file, output_file, lines_to_read)
