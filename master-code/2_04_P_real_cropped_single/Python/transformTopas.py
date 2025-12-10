import pandas as pd
import numpy as np

def convert_voxel_to_trajectory(input_file, output_file, filtered_output_file, max_y=224):
    # Read the CSV file, skipping unnecessary header lines
    df = pd.read_csv(input_file, comment='#', header=None)
    
    # Extract the columns
    df.columns = ['x', 'y', 'z', 'dose']
    
    # Filter out rows with dose value 0
    df_filtered = df[df['dose'] > 0].copy()
    
    # Calculate the actual y position by mirroring the y values using .loc to avoid SettingWithCopyWarning
    df_filtered.loc[:, 'y_actual'] = max_y - df_filtered['y']
    
    # Select the required columns and rename them
    df_trajectory = df_filtered[['x', 'y_actual', 'dose']]
    df_trajectory.columns = ['x', 'y', 'dose']
    
    # Save the new trajectory to a CSV file without a header
    df_trajectory.to_csv(output_file, index=False, header=False)
    
    # Filter further for the dose greater than 1e-5
    df_filtered_dose = df_trajectory[df_trajectory['dose'] > 1e-5]
    
    # Save the filtered trajectory to a new CSV file without a header
    df_filtered_dose.to_csv(filtered_output_file, index=False, header=False)

# Convert the TOPAS_Beam.csv to trajectory_TOPAS.csv and trajectory_TOPAS_filtered.csv
convert_voxel_to_trajectory('../TOPAS/TOPAS_Beam.csv', 'trajectory_TOPAS.csv', 'trajectory_TOPAS_filtered.csv')
#convert_voxel_to_trajectory('TOPAS_Beam_Slice1.csv', 'trajectory_TOPAS.csv', 'trajectory_TOPAS_filtered.csv')

