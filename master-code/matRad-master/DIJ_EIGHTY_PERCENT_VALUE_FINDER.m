%% EACH BIXEL DOSE DISTRIBUTION

% Dimensions of the dose grid (replace these with actual values)
x_dim = dij.doseGrid.dimensions(2); 
y_dim = dij.doseGrid.dimensions(1);  
z_dim = dij.doseGrid.dimensions(3);  

% Reshape each column of physicalDoseData into a 3D grid
for i = 1:size(dij.physicalDose{1,1}, 2)  % Iterate through each bixel
    % Convert sparse column to full matrix
    fullDoseData = full(dij.physicalDose{1,1}(:, i));  % Convert the sparse data to full
    
    % Reshape the full dose data into a 3D grid
    dose_3D = reshape(fullDoseData, [y_dim, x_dim, z_dim]);
    
    % Extract the 2D slice corresponding to z=2
    z_slice = 2;  % Slice number for z = 2
    slice_2_data = dose_3D(:, :, z_slice);
    
    % Save the 2D slice to a CSV file (one per bixel)
    csv_filename = sprintf('physicalDoseData_slice_z_%d_bixel_%d.csv', z_slice, i);
    writematrix(slice_2_data, csv_filename);
    
    % Confirmation message
    disp(['2D dose distribution for z=', num2str(z_slice), ' for bixel ', num2str(i), ' has been saved to ', csv_filename]);
end


%% MAX DOSE VALUES

% Initialize maxRows to store the 1D indices of maximum doses
maxRows = zeros(dij.totalNumOfBixels, 1);

% Initialize maxRows_3D to store the 3D coordinates corresponding to each 1D index
maxRows_3D = zeros(dij.totalNumOfBixels, 3);  % Each row will hold [x, y, z]

% Dimensions of the 3D matrix
x_dim = dij.doseGrid.dimensions(2);  % Number of columns
y_dim = dij.doseGrid.dimensions(1);  % Number of rows
z_dim = dij.doseGrid.dimensions(3);  % Number of slices

% Number of elements per slice (x * y)
elements_per_slice = x_dim * y_dim;

% Track skipped entries
skipped_entries = [];
valid_indices = [];

for i = 1:dij.totalNumOfBixels
    % Find the 1D index of the maximum dose for the current bixel
    maxRowsRaw = find(dij.physicalDose{1, 1}(:, i) == max(dij.physicalDose{1, 1}(:, i)));
    maxRows(i) = maxRowsRaw;  % Store the 1D index

    % Convert the 1D index to 3D coordinates (x, y, z)
    linear_index = maxRowsRaw;
    
    % Determine z index
    z = floor((linear_index - 1) / elements_per_slice) + 1;
    
    % Calculate remaining index within the slice
    remainder = linear_index - (z - 1) * elements_per_slice;
    
    % Determine x index
    x = floor((remainder - 1) / y_dim) + 1;
    
    % Determine y index
    y = remainder - (x - 1) * y_dim;

    % Store the 3D coordinates in maxRows_3D
    maxRows_3D(i, :) = [x, y, z];
    
    % If z is not equal to 2, track the skipped entry
    if z ~= 2
        skipped_entries = [skipped_entries; i, x, y, z];  % Track index and coordinates
    else
        % Track valid indices where z = 2
        valid_indices = [valid_indices; i];
    end
end

% Filter entries where z equals 2
maxRows_3D_secondSlice = maxRows_3D(maxRows_3D(:, 3) == 2, :);

% Save the filtered entries to a CSV file
writematrix(maxRows_3D_secondSlice, 'maxRows3D_secondSlice.csv');

% Print skipped entries with z other than 2
if ~isempty(skipped_entries)
    disp('Skipped entries where z is not equal to 2:');
    disp(array2table(skipped_entries, 'VariableNames', {'EntryIndex', 'X', 'Y', 'Z'}));
else
    disp('No skipped entries, all z values are equal to 2.');
end

% Confirmation message
disp('Filtered entries with z = 2 have been saved to maxRows3D_secondSlice.csv');


%% ENERGY

% Initialize an empty array to store all energy values
all_energy_values_1D = [];

% Loop through each ray (1x8 struct)
for i = 1:length(stf.ray)
    % Access the energy values for the current ray
    energy_values = stf.ray(i).energy;
    
    % Concatenate the energy values into the 1D array
    all_energy_values_1D = [all_energy_values_1D; energy_values(:)];
end

% Save all energy values to a CSV file in 1D format (one column)
writematrix(all_energy_values_1D, 'all_energy_values_1D.csv');

% Confirmation message
disp('All energy values have been saved to all_energy_values_1D.csv in 1D format.');

% Remove the skipped entries (where z ~= 2) from energy values
energy_values_secondSlice = all_energy_values_1D(valid_indices);

% Save the filtered energy values to a CSV file in 1D format (one column)
writematrix(energy_values_secondSlice, 'energy_values_secondSlice.csv');

% Confirmation message
disp('Filtered energy values for the second slice (z = 2) have been saved to energy_values_secondSlice.csv.');


%% ENERGY AND POSITION

% Load or use the existing maxRows_3D_secondSlice and energy_values_secondSlice

% Ensure maxRows_3D_secondSlice and energy_values_secondSlice have compatible sizes
if size(maxRows_3D_secondSlice, 1) == length(energy_values_secondSlice)
    % Replace the third column of maxRows_3D_secondSlice with the energy values
    maxRows_3D_with_energy = maxRows_3D_secondSlice;  % Make a copy
    maxRows_3D_with_energy(:, 3) = energy_values_secondSlice;  % Replace third column with energy
    
    % Save the updated maxRows_3D_with_energy to a CSV file
    writematrix(maxRows_3D_with_energy, 'maxRows_3D_with_energy.csv');
    
    % Confirmation message
    disp('The third column of maxRows_3D_secondSlice has been replaced with energy values and saved to maxRows_3D_with_energy.csv');
else
    error('Size mismatch: The number of rows in maxRows_3D_secondSlice and the number of energy values in energy_values_secondSlice must be the same.');
end

%% DOSE ALONG X-AXIS (DEPTH DOSE CURVE)

% Assuming dose_3D is already computed and we are using z_slice = 2.
% Sum dose along the y-dimension for each x across all y values in the slice.

% Initialize a vector to store the summed dose values along x
dose_along_x = sum(slice_2_data, 1);  % Summing over the y-dimension

% Plot the depth-dose curve along the x-axis
figure;
plot(1:x_dim, dose_along_x, '-o');  % x-axis goes from 1 to x_dim, plot dose_along_x
title('Depth-Dose Curve along x-axis');
xlabel('x (depth)');
ylabel('Dose');
grid on;

