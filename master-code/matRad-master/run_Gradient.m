% Load data and initialize variables as before
data = readmatrix('x_y_energy_values_at_80_percent.csv');
y_positions = data(:, 2);  % Column 2 represents y-position
energy_values = data(:, 3);  % Column 3 represents energy values
new_x_positions = zeros(size(y_positions));
new_y_positions = zeros(size(y_positions));  % To store updated y-positions if adjusted
new_energy_values = zeros(size(energy_values));  % To store updated energy values

% Initialize arrays to store initial positions and adjusted starting positions
adjusted_starting_x = zeros(size(y_positions));
adjusted_starting_y = zeros(size(y_positions));
adjusted_starting_energy = zeros(size(energy_values));

% Initialize arrays to store the initial positions under the magnetic field (for verification)
initial_x_positions_B = zeros(size(y_positions));
initial_y_positions_B = zeros(size(y_positions));
initial_energy_values_B = zeros(size(energy_values));

magnetic_field = [0; 0; 1.5]; % Magnetic field in Tesla
SPR = readmatrix('SPR_values.csv');

% Parameters for gradient descent
learning_rate_E = 0.0005; % Step size for energy
learning_rate_y = 0.0005; % Step size for position (if y is adjusted)
tolerance = 1e-4;         % Convergence tolerance
max_steps = 200;          % Maximum iterations

% Error threshold definitions
small_error_threshold = 0.5;  % Reduce small difference threshold for finer control
big_error_threshold = 1.5;    % Reduce big difference threshold for finer control

% Loop over each row in the CSV data
for i = 1:length(y_positions)
    % Current y and energy
    y_position = y_positions(i) / 10.9375;  % Divide y-position by pixel factor
    E0_MeV = energy_values(i);  % Energy value

    % Target without magnetic field
    Target = analyzeProtonTrajectory(y_position, E0_MeV, [0; 0; 0], SPR);

    % Initial position with magnetic field (Target_start)
    Target_start = analyzeProtonTrajectory(y_position, E0_MeV, magnetic_field, SPR);

    % Save the initial positions under magnetic field for verification (before gradient descent)
    initial_x_positions_B(i) = Target_start(1) * 10;  % Scale back to original x-position
    initial_y_positions_B(i) = Target_start(2) * 10.9375;  % Scale back to original y-position
    initial_energy_values_B(i) = E0_MeV;  % Initial energy before gradient descent

    % Gradient Descent Loop
    Ei = E0_MeV;  % This energy will be adjusted
    yi = y_position;  % This y-position will be adjusted

    % Initialize variables to track the starting positions during the gradient descent
    initial_yi = yi;  % Keep track of the starting y_position
    initial_Ei = Ei;  % Keep track of the starting energy

    for step_i = 1:max_steps
        % Compute the difference between the current position and target
        position_this_step = analyzeProtonTrajectory(yi, Ei, magnetic_field, SPR);
        error_x = position_this_step(1) - Target(1);
        error_y = position_this_step(2) - Target(2);

        % Check for convergence
        if abs(error_x) < tolerance && abs(error_y) < tolerance
            break;
        end

        % Adjust learning rates based on error size
        if abs(error_x) >= big_error_threshold
            % Larger errors require larger steps, but within smaller thresholds now
            dE = -error_x * E0_MeV * learning_rate_E;
            dy = error_x * y_position * learning_rate_y;
        elseif abs(error_x) <= small_error_threshold
            % For small errors, fine-tune the adjustments
            dE = -error_x * E0_MeV * learning_rate_E * 0.5;
            dy = error_x * y_position * learning_rate_y * 0.5;
        end

        % Update energy and position using gradient descent
        Ei = Ei + dE; % Update energy
        yi = yi + dy; % Update y-position
    end

    % Save the final x-position, updated y-position, and energy after gradient descent
    new_x_positions(i) = position_this_step(1) * 10;  % Convert x to original scale
    new_y_positions(i) = yi * 10.9375;  % Convert back to original scale for saving
    new_energy_values(i) = Ei;  % Store the updated energy

    % Calculate the difference in y-position, x-position, and energy between Target_start and final values
    adjusted_starting_x(i) = Target_start(1) - position_this_step(1);  % Difference in x
    adjusted_starting_y(i) = Target_start(2) - yi;  % Difference in y-position
    adjusted_starting_energy(i) = E0_MeV - Ei;  % Difference in energy
end

% Combine new x-positions, updated y-positions, and updated energy values into a matrix
new_data = [new_x_positions, new_y_positions, new_energy_values];
writematrix(new_data, 'x_y_energy_values_at_80_percent_gradient.csv');

% Save the initial starting positions before gradient descent under the magnetic field (Target_start) in a CSV file
initial_B_data = [initial_x_positions_B, initial_y_positions_B, initial_energy_values_B];
writematrix(initial_B_data, 'x_y_energy_values_at_80_percent_gradient_B.csv');

% Save the adjusted starting positions (differences in x, y, and energy to reach the target) in a CSV file
adjusted_starting_data = [adjusted_starting_x * 10, adjusted_starting_y * 10.9375, adjusted_starting_energy];
writematrix(adjusted_starting_data, 'x_y_new_starting.csv');


%%%%%%%%%%%%%%%%%% Function to Analyze Proton Trajectory %%%%%%%%%%%%%%%%%%

function final_position = analyzeProtonTrajectory(y, E, B, SPR)
    % Constants (needed for velocity calculation for example)
    atomic_mass_unit_MeV_c2 = 931.494; % MeV/c^2
    c = 299792458 * 100; % [mm/s]
    
    grid_step = 0.00109375; % [m] This is the CT grid step /"PixelSpacing" 
    total_distance = 83.91796875; % [mm] CT rows * CT columns * CT grid step
    
    % Pre-calculations for initial velocity
    gamma0 = 1 + E / atomic_mass_unit_MeV_c2;
    beta0 = sqrt(1 - 1/gamma0^2);
    v0 = beta0 * c;
    
    % Set initial conditions based on input parameters
    initial_position = [0; y; 0]; % Initial y-position
    initial_velocity = [v0; 0; 0]; % Initial velocity
    
    % Create an instance of ProtonSimulation
    protonSim = ProtonSimulation_half(SPR, E, B, grid_step, total_distance, initial_position, initial_velocity);
    
    % Initialize and run the simulation
    protonSim = protonSim.initializeStep();
    protonSim = protonSim.simulate();
    protonSim.saveResults();

    % Get the final position using the number of non-zero steps
    num_non_zero_steps = protonSim.displayStepRange();
    if num_non_zero_steps > 0
        final_position = protonSim.positions(:, num_non_zero_steps); % Extract the last valid position
    else
        final_position = [0; 0; 0]; % Default if no valid positions found
    end
end
