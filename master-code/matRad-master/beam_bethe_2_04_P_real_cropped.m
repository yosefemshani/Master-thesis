%% Integrating magnetic field grid size model with consideration of proton energy loss applied on a given CT patient data set (2_04_P_real_cropped)
%% Author: Yosef Emshani, 01.03.2024 17:51


%% Clear workspace

clear;  % Clear workspace
clc;    % Clear command window
close all;  % Close all figures

%% Read interpolation data from https://libamtrack.github.io/web_dev/#Particlerangeandstoppingpower/Stoppingpowervsenergy
% Used start energy 0.025 MeV, end Energy 1000 MeV, 500 Points, Liquid water, proton, ICRU

dEdx_data = csvread('energy_interpolation_ICRU.csv');

energy_data = dEdx_data(:,1);
stopping_power_data = dEdx_data(:,2);


%% Definition of constant variables

mp_MeV = 938.272013;
mp = 1.67262158e-27; % [kg]
c = 299792458*100; % [mm/s]
atomic_mass_unit_MeV_c2 = 931.494; %MeV/c^2; double check this value later once code runs! might be 931.494e6... 
q = 1.602e-19; % elementary charge: positive value for protons!


%% Definition of variables one might change

E0_MeV = 100; % initial proton energy, 332 MeV drops to 0, 333 doesnt..
B = [0; 0; 1.5]; % [T]
desired_final_position = [10; 15; 0]; % [cm]


%% Definition of Lorentz force in 2D
% Let us say that v = (vx, vy, 0 [since we only want 2 dimensions]) and
% B = (0, 0, B0 [with any B0, we couldd use 1.5 T for example]);
% here we also asumme that B is only in z-direction!
% Since the Lorentz force is F = q * (v x B), we would calculate v x B:
% v x B = (vx,vy,0) x (0,0,B0) = (0, 0, vxB0) - (0, 0, vyB0)
% = (-vyB0, vxB0, 0)
% Thus in 2D, we receive: F = q * (-vyB0, vxB0, 0) or
% Fx = -q * vy * B0 and Fy = q * vx * B0

%% Initialize first step of velocity and acceleration with Lorentz force
% relativistic consideration!

gamma0 = 1 + E0_MeV / atomic_mass_unit_MeV_c2;  % lorentz factor for initial energy
beta0 = sqrt(1 - 1/gamma0^2); % initial relativistic speed beta = v/c
v0 = beta0 * c; % initial velocity relative to the speed of light
disp(['initial velocity is : ', num2str(v0)]);

v = [v0; 0; 0]; % Velocity vector with velocity in x direction!

%Fx = -q * vy * B0;
%Fy = q * vx * B0;

F = q * cross(v, B);

%ax = Fx/mp_MeV; 
%ay = Fy/mp_MeV;

a = F/mp;


%% Considerations for simulation
% This CT data set has a pixel spacing (vertical and horizonal) of 1.09375 [mm]
% This leads to a grid step size of 0.00109375 m!
% Also, there are 225 total rows and 341 total columns
% Finally, there is a total of 76725 pixels!
% Do not forget that these values are for ONE slice, not all slices!!

% Total distance horizonally is total pixels * pixel spacing
% Thus, distance: 76725 * 1.09375 mm = 83917.96875 mm = 83.91796875 m
% seems like a large distance for a patient... DOUBLE CHECK! 

% Rescaling the grid size from given pixel spacing 1.09375 [mm] to 0.001 mm
% This needs to be done by calculating the new total distance
% new total distance = total distance * (new grid step / old grid step)
% new total distance = 83.91796875 m * (0.00001 m / 0.00109375 m)
% new total distance ≈ 0.76725 m
% now we can use 0.00001 as a grid step size!

% With rescaled grid size of 0.01 mm pixel spacing, we get
% total distance 76725 * 0.001 mm = 76.725 mm = 0.076725 m

%% Initialization for simulation
%grid_step = 0.001; % [mm]
%total_distance = 76.725; % [mm]
grid_step = 0.00109375; % [m]
total_distance = 83.91796875; % [mm]
steps = round(total_distance / grid_step); % number of steps to update

% Define number of steps
num_steps = steps;

% Initialize arrays to store our results
positions = zeros(3, num_steps); % 2D array to store position vectors
velocities = zeros(3, num_steps); % 2D array to store velocity vectors
accelerations = zeros(3, num_steps); % 2D array to store acceleration vectors
times = zeros(3, num_steps); % 2D array to store times needed for x-y position to get to grid step

% Initial conditions
positions(:,1) = [0; 11.25; 0]; % Initial position (assuming starting from origin)
%velocities(:,1) = [vx; vy]; % Initial velocity
velocities(:,1) = v;

% Initialize proton's energy
E_MeV = E0_MeV;

% Initialize x_distance
x_distance = 0;

% Import the SPR values from the CSV file
SPR_values = readmatrix('water.csv');%readmatrix('IMG0035_cropped_SPR.csv');

% Ensure SPR_values has the correct dimensions
assert(numel(SPR_values) == num_steps, 'Number of SPR values does not match the number of steps');

% Initialize array to store energy losses adjusted by SPR
stopping_power_adjusted = zeros(1, num_steps);

%% Simulation
% Loop through each step
for step = 2:num_steps
    % Calculate time taken to travel each step
    time_step = grid_step / norm(velocities(:, step - 1));

    F = q * cross((velocities(:, step - 1)), B);
    a = F/mp;

    % Check if energy is below the BETHE_LOWER_LIMIT
    if E_MeV < 0.025
        % Set energy loss to 0
        stopping_power_interp = 0;
        % Set energy to 0
        E_MeV = 0;
        % Break out of the loop
        disp(['Energy dropped below Bethe lower limit, E=0 MeV']);
        break;
    else
        % Calculate energy loss without SPR adjustment
        stopping_power_interp = interp1(energy_data, stopping_power_data, E_MeV, 'linear'); % [keV/um]
        stopping_power_interp = stopping_power_interp * 10^(1); % [MeV/mm] think about this later

        % Adjust energy loss by multiplying with SPR value for this step
        stopping_power_interp_adjusted = stopping_power_interp * SPR_values(step);
        
        % Store adjusted energy loss for this step
        stopping_power_adjusted(step) = stopping_power_interp_adjusted;
        
        % Update energy for the next step
        E_MeV = E_MeV - stopping_power_interp_adjusted * grid_step;
    end

    % Recalculate gamma and beta
    gamma = 1 + E_MeV / atomic_mass_unit_MeV_c2;
    % Ensure gamma is scalar
    gamma = gamma(1);
    beta = sqrt(1 - 1/gamma^2);

    % Check if gamma drops below 1
    if gamma < 1
        disp(['Gamma dropped below 1 at step ', num2str(step)]);
        break;
    end
    
 % Update velocity
    vx_new = velocities(1, step-1) + a(1) * time_step;
    vy_new = velocities(2, step-1) + a(2) * time_step;
    vz_new = velocities(3, step-1) + a(3) * time_step;
    velocities(:, step) = [vx_new; vy_new; vz_new];

% Update position based on distance traveled in x-direction
x_distance = x_distance + abs(vx_new) * time_step;
if x_distance >= grid_step
    % Calculate the number of full steps in x-direction
    full_steps = floor(x_distance / grid_step);
    
    % Update the x component based on full steps
    x_new = positions(1, step - 1) + sign(vx_new) * full_steps * grid_step;
    
    % Update remaining distance after full steps
    x_distance = x_distance - full_steps * grid_step;
    
    % Update y and z components accordingly
    y_new = positions(2, step - 1) + vy_new * (full_steps * grid_step / abs(vx_new));
    z_new = positions(3, step - 1) + vz_new * (full_steps * grid_step / abs(vx_new));
    
    % Update positions
    positions(:, step) = [x_new; y_new; z_new];
else
    % If the distance traveled is less than grid_step, update only x component
    x_new = positions(1, step - 1) + sign(vx_new) * grid_step;
    % x_new = positions(1, step - 1) + sign(velocities(1, step - 1)) * full_steps * grid_step;
    y_new = positions(2, step - 1) + vy_new * (x_distance / abs(vx_new));
    z_new = positions(3, step - 1) + vz_new * (x_distance / abs(vx_new));
    
    % Update positions
    positions(:, step) = [x_new; y_new; z_new];
    
    % Update remaining distance after the step
    x_distance = 0;
end

    % Store time needed for x-y position to get to grid step
    times(step) = times(step - 1) + time_step;
    
    % Store accelerations
    accelerations(:, step) = a;

    disp(['dE/dx is ', num2str(stopping_power_interp), ' [MeV/mm]']);
    disp(['current Energy is ', num2str(E_MeV), ' [MeV]']);
    disp(['---------------']);

end

% Find the index of the last non-zero position
last_nonzero_index = find(sum(positions ~= 0, 1), 1, 'last');

% Get the last non-zero position
last_position = positions(:, last_nonzero_index);

% Calculate the shift required
shift = desired_final_position - last_position;

% Apply shift to all positions
positions = positions + shift;

% Check if final position has been reached
if isequal(positions(:, last_nonzero_index), desired_final_position)
    % Set all subsequent positions to zero
    positions(:, last_nonzero_index+1:end) = 0;
end

% Save positions to CSV file after the loop completes
csvwrite('trajectory.csv', positions');
disp(['Positions have been saved!']);


% Plot trajectory
figure;
plot(positions(1, 1:step-1), positions(2, 1:step-1), 'b.-', 'LineWidth', 1);
title('Trajectory of Proton under Magnetic Field with Energy Loss');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
grid on;
axis equal;
