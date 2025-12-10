% run_Loop.m

% Load the original CSV file with x, y, and energy values
data = readmatrix('x_y_energy_values_at_80_percent.csv');

% Load the new CSV file with x, y, and energy shift values
data_shift = readmatrix('x_y_new_starting.csv');

% Extract y-positions and energy values from columns 2 and 3 of the original data
y_positions = data(:, 2);  % Column 2 represents y-position
energy_values = data(:, 3);  % Column 3 represents energy values

% Initialize array to store new x-positions and y-positions
new_x_positions = zeros(size(y_positions));
new_y_positions = zeros(size(y_positions));

% Magnetic field (kept constant)
magnetic_field = [0; 0; 1.5]; % Example magnetic field in T

SPR = readmatrix('SPR_values.csv'); 
% Load SPR values SPR_values(_slice).csv watertestSPR.csv bone_SPR.csv heavybone_SPR.csv

% Loop over each row in the CSV data
for i = 1:length(y_positions)
    % Current y and energy
    y_position = y_positions(i) / 10.9375;  % Divide y-position by 10 for calculation and 1.09375 for pixel
    E0_MeV = energy_values(i);  % Energy value

    % Adjust initial conditions using values from the data_shift CSV
    % new initial position for this step: [0 - x_shift, y - y_shift, E - E_shift]
    initial_position = [0 - data_shift(i, 1); y_position - data_shift(i, 2); E0_MeV - data_shift(i, 3)];

    % Calculate the final x-position using analyzeProtonTrajectory function
    final_pos = analyzeProtonTrajectory(initial_position, E0_MeV, magnetic_field, SPR);

    % Multiply the resulting x-position by 10 before saving
    new_x_positions(i) = final_pos(1) * 10;

    % Multiply the resulting y-position by 10.9375 before saving
    new_y_positions(i) = final_pos(2) * 10.9375;
end

% Combine new x-positions, adjusted y-positions, and energy values into a matrix
new_data = [new_x_positions, new_y_positions, energy_values];

% Save the new matrix into a new CSV file
writematrix(new_data, 'x_y_new_MATLAB.csv');

% Updated analyzeProtonTrajectory function to accept the adjusted initial position
function final_position = analyzeProtonTrajectory(initial_position, E, B, SPR)
    % Constants (needed for velocity calculation for example)
    atomic_mass_unit_MeV_c2 = 931.494; % MeV/c^2
    c = 299792458 * 100; % [mm/s]
    
    grid_step = 0.00109375; % [m] This is the CT grid step /"PixelSpacing"
    total_distance = 83.91796875; % [mm] CT rows * CT columns * CT grid step
    
    % Pre-calculations for initial velocity
    gamma0 = 1 + E / atomic_mass_unit_MeV_c2;
    beta0 = sqrt(1 - 1/gamma0^2);
    v0 = beta0 * c;
    
    % Initial velocity stays the same as in the original code
    initial_velocity = [v0; 0; 0]; % Initial velocity
    
    % Create an instance of ProtonSimulation
    protonSim = ProtonSimulation_half(SPR, E, B, grid_step, total_distance, initial_position, initial_velocity);
    
    % Initialize and run the simulation
    protonSim = protonSim.initializeStep();
    protonSim = protonSim.simulate();
    protonSim.saveResults();

    % Display the number of steps with non-zero values
    % Get the number of steps with non-zero values
    num_non_zero_steps = protonSim.displayStepRange();
            
    % Get the final position using the number of non-zero steps
    if num_non_zero_steps > 0
       final_position = protonSim.positions(:, num_non_zero_steps); % Extract the last valid position
       %disp(['Final valid position (x, y, z): ', num2str(final_position')]);
    else
       %disp('No valid positions found.');
    end
end
