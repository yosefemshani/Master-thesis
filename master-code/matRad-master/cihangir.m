% Initialization for storing x, y values, 80% dose values, and energy
x_y_energy_values_at_80_percent = [];

% Counter for the current bixel (across all rays)
current_bixel_id = 1;

% Loop through all rays
for ray_id = 1:length(stf.ray)
    % Get the x and y position of the current ray relative to the isocenter
    ray_y_pos = stf.ray(ray_id).rayPos(2);
    
    % Get the isocenter
    isoCenter_y = stf.isoCenter(2);

    % Actual y position of the ray (relative to the isocenter)
    ray_y_absolute = isoCenter_y + ray_y_pos;

    % Number of bixels for this ray
    num_bixels = stf.numOfBixelsPerRay(ray_id);
    
    % Loop through all bixels in this ray
    for bixel_idx = 1:num_bixels
        % Step 1: Extract the 1D vector for the current bixel
        doseVector = full(dij.physicalDose{1,1}(:, current_bixel_id));
        LETVector = full(dij.mLETDose{1,1}(:, current_bixel_id));

        % Dimensions of the dose grid
        x_dim = dij.doseGrid.dimensions(2); % Number of x coordinates
        y_dim = dij.doseGrid.dimensions(1); % Number of y coordinates
        z_dim = dij.doseGrid.dimensions(3); % Number of z layers

        % Step 2: Convert the 1D vector into a 3D matrix
        doseMatrix = reshape(doseVector, [y_dim, x_dim, z_dim]);
        LETMatrix = reshape(LETVector, [y_dim, x_dim, z_dim]);

        % Step 3: Restrict to the second z layer
        doseMatrix_z2 = doseMatrix(:,:,2); % Dose values in the second z layer
        LETMatrix_z2 = LETMatrix(:,:,2); % LET values in the second z layer

        % Step 4: Sum the dose values along the y-axis for each x coordinate
        depthDoseCurve = sum(doseMatrix_z2, 1); % Sum over y for each x coordinate
        depthLETCurve = sum(LETMatrix_z2, 1); % Sum over y for each x coordinate

        % Step 5: Normalize the dose values to the maximum value (relative dose)
        max_dose = max(depthDoseCurve);
        relativeDepthDoseCurve = depthDoseCurve / max_dose;

        max_LET = max(depthLETCurve);
        relativeDepthLETCurve = depthLETCurve / max_LET;

        % Step 6: Find the x value where the dose falls to 80% of the maximum (after the peak)
        dose_threshold = 0.8; % 80% of the maximum value
        % Find the peak
        [~, max_index] = max(relativeDepthDoseCurve);
        % Find the first value after the peak that falls below 80%
        idx_above_threshold = find(relativeDepthDoseCurve(max_index:end) >= dose_threshold) + max_index - 1;
        idx_below_threshold = find(relativeDepthDoseCurve(max_index:end) < dose_threshold, 1, 'first') + max_index - 1;

        % Interpolation between the two points
        if ~isempty(idx_above_threshold) && ~isempty(idx_below_threshold) && idx_below_threshold > idx_above_threshold(end)
            x1 = idx_above_threshold(end); % Last index where dose >= 80%
            x2 = idx_below_threshold; % First index where dose < 80%
            D1 = relativeDepthDoseCurve(x1); % Dose at x1
            D2 = relativeDepthDoseCurve(x2); % Dose at x2

            % Linear interpolation for the exact x value where the dose falls to 80%
            x_80_percent = x1 + (dose_threshold - D1) / (D2 - D1) * (x2 - x1);
        else
            x_80_percent = NaN; % If no suitable point is found
        end
        
        % figure;
        % hold on;
        % plot((1:x_dim)*1.09375, depthLETCurve/max_LET, 'b.-', 'LineWidth', 1);
        % plot((1:x_dim)*1.09375, depthDoseCurve/max_dose, 'r.-', 'LineWidth', 1);
        % hold off;
        % xlabel('Range (mm)');
        % ylabel('mLET vs Dose');
        % grid on;
        % ylim([0,1]);

        % Extract the energy of the current bixel from the current ray
        energy = stf.ray(ray_id).energy(bixel_idx);
        
        % Store the x value (at 80% dose), the corrected y value (Ray + Isocenter), and the energy
        x_y_energy_values_at_80_percent = [x_y_energy_values_at_80_percent; x_80_percent*dij.ctGrid.resolution.x, ray_y_absolute, energy];
        
        % Increment the bixel counter
        current_bixel_id = current_bixel_id + 1;
    end
end

% Define the filename
%filename = 'C:/Users/MP_Bernoulli/Downloads/master-code/matRad-master/x_y_energy_values_at_80_percent.csv';
filename = 'x_y_energy_values_at_80_percent.csv';

% Save the matrix to a CSV file
writematrix(x_y_energy_values_at_80_percent, filename);