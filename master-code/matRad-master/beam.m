%% Trajectory of a single proton with specified intial energy value under influence of a magnetic field with consideration of a grid size
%% Author: Yosef Emshani, 23.02.2024 13:40

%% Clear workspace

clear;  % Clear workspace
clc;    % Clear command window
close all;  % Close all figures


%% Definition of constant variables

mp_MeV = 938.272013;
c = 299792458; % [m/s]
atomic_mass_unit_MeV_c2 = 931.494; %MeV/c^2; double check this value later once code runs! might be 931.494e6... 
q = 1.602e-19; % elementary charge: positive value for protons!


%% Definition of variables one might change

E0_MeV = 70; % initial proton energy
B0 = 1.5; % [T]


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

vx = v0; % assuming initial velocity is only headed towards x-direction
vy = 0.2*c; 

Fx = -q * vy * B0;
Fy = q * vx * B0;

ax = Fx/mp_MeV; 
ay = Fy/mp_MeV;

%% Considerations for simulation
% Let us assume a 10x10cm^2 grid with grid step size 1cm
% Let us also assume that energy, velocity, force and acceleration are only
% updated once that proton travels 1cm!
% We will need to calculate the time, since (instead of given step time, we
% have given step size) and with that time we can calculate our variables!

grid_step = 0.01; % [m]
total_distance = 1; % [m]
steps = round(total_distance / grid_step); % number of steps to update

% Define number of steps
num_steps = steps;

% Initialize arrays to store our results
positions = zeros(2, num_steps); % 2D array to store position vectors
velocities = zeros(2, num_steps); % 2D array to store velocity vectors
accelerations = zeros(2, num_steps); % 2D array to store acceleration vectors
times = zeros(2, num_steps); % 2D array to store times needed for x-y position to get to grid step

% Initial conditions
positions(:,1) = [0; 0]; % Initial position (assuming starting from origin)
velocities(:,1) = [vx; vy]; % Initial velocity

% Define number of steps
num_steps = steps;

% Loop through each step
x_distance = 0;  % Initialize distance traveled in x-direction
for step = 2:num_steps
    % Calculate time taken to travel each step
    time_step = grid_step / norm(velocities(:, step - 1));
    
    % Update velocity components using Lorentz force equation
    Fx = -q * velocities(2, step - 1) * B0;
    Fy = q * velocities(1, step - 1) * B0;
    ax = Fx / mp_MeV;
    ay = Fy / mp_MeV;
    vx_new = velocities(1, step - 1) + ax * time_step;
    vy_new = velocities(2, step - 1) + ay * time_step;
    velocities(:, step) = [vx_new; vy_new];
    
    % Update distance traveled in x-direction
    x_distance = x_distance + abs(vx_new) * time_step;
    
    % Update position based on distance traveled in x-direction
    if x_distance >= grid_step
        % Update position only if x-coordinate has moved by 1 cm
        x_new = positions(1, step - 1) + sign(vx_new) * grid_step;
        y_new = positions(2, step - 1) + vy_new * (grid_step / abs(vx_new));
        positions(:, step) = [x_new; y_new];
        x_distance = x_distance - grid_step;  % Reset distance traveled in x-direction
    else
        % Keep previous position
        positions(:, step) = positions(:, step - 1);
    end
    
    % Store time needed for x-y position to get to grid step
    times(step) = times(step - 1) + time_step;
    
    % Store accelerations
    accelerations(:, step) = [ax; ay];
end


% Plot trajectory
figure;
plot(positions(1,:), positions(2,:), 'b.-', 'LineWidth', 1);
title('Trajectory of Proton under Magnetic Field');
xlabel('X Position (m)');
ylabel('Y Position (m)');
grid on;
axis equal;
