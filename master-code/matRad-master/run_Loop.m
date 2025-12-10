% run_Loop.m

% Load the original CSV file with x, y, and energy values
data = readmatrix('x_y_energy_values_at_80_percent.csv');

% Extract y-positions and energy values from columns 2 and 3
y_positions = data(:, 2);  % Column 2 represents y-position
energy_values = data(:, 3);  % Column 3 represents energy values

% Initialize array to store new x-positions
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

    % Calculate the final x-position using analyzeProtonTrajectory function
    final_pos = analyzeProtonTrajectory(y_position, E0_MeV, magnetic_field, SPR);

    % Multiply the resulting x-position by 10 before saving
    new_x_positions(i) = final_pos(1) * 10;

    % Multiply the resulting y-position by 10 before saving
    new_y_positions(i) = final_pos(2) * 10.9375;
end

% Multiply y-positions by 1.09375 (for saving, not for calculations)
% adjusted_y_positions = y_positions * 1.09375;

% Combine new x-positions, adjusted y-positions, and energy values into a matrix
new_data = [new_x_positions, new_y_positions, energy_values];

% Save the new matrix into a new CSV file
writematrix(new_data, 'x_y_energy_values_at_80_percent_MATLAB_shifted.csv');





%%%%%%%%%%%%%%%%%% Liheng write here
%
%
%Target = analyzeProtonTrajectory(y_position, E0_MeV, [0;0;0], SPR);
%Target_start = analyzeProtonTrajectory(y_position, E0_MeV, magnetic_field, SPR);
%
%Ei = E0_MeV;
%yi = y_position;
%
%%%%% those are constant
%step_size_E=0.01;
%step_size_y=0.01;
%lambda_E = E0_MeV*(Target_start(1)-Target(1))/Target(1);
%lambda_y = y_position*(Target_start(1)-Target(1))/Target(1);
%
%%%%% change of the first step
%dE= (Target_start(1)-Target(1))*lambda_E;
%dy= -(Target_start(1)-Target(1))*lambda_y;
%
%for stpe_i = 1:100
%    Ei = Ei+dE;
%    yi=yi+dy;
%    position_this_step = analyzeProtonTrajectory(yi, Ei, magnetic_field, SPR);
%    %disp((Target(1)-position_this_step(1))^2 + (Target(2)-position_this_step(2))^2)
%    %disp(Target(1)-position_this_step(1))
%    dE= (position_this_step(1)-Target(1))*lambda_E;
%    dy= -(position_this_step(1)-Target(1))*lambda_y;
%  %  break
%
%
%end



%%%%%%%%%%%%%%%%%% Liheng write here



% New Function to Analyze Proton Trajectory
function final_position = analyzeProtonTrajectory(y, E, B, SPR)
    % Constants (needed for velocity calculation for example)
    atomic_mass_unit_MeV_c2 = 931.494; % MeV/c^2
    c = 299792458 * 100; % [mm/s]
    
    grid_step = 0.00109375; % [m] This is the CT grid step /"PixelSpacing"  one can aquire by looking at metadata using 3D Slicer for example
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
