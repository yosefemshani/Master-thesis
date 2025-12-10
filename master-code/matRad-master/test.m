% analysis of Lihengs idea on reverse dose calculation method with given
% end position, velocity and angle

% definition of basic 3D 10x10x10 cm^3 box
x = 10;
y = 10;
z = 10;

% definition of parameters
q = 1.602 * 10^(-19); % [C]
m = 1.6726 * 10^(-27); % [kg]
B = [0, 1.5, 0]; % 1.5 [T] in y-axis

% initial conditions
initial_position = [0,0,0]; % initial position of x y and z in [m]
initial_velocity = [269813212.2, 0, 0]; % initial velocity v = 0.9c in x-axis [m/s]


% simulation parameters
dt = 1 * 10^(-9); % time step in [s]
num_steps = 1000; % 1000 steps in iteration

% Initialize arrays to store position, velocity, and angles at each step
positions = zeros(num_steps + 1, 3);
velocities = zeros(num_steps + 1, 3);
angles = zeros(num_steps + 1, 2); % Columns represent azimuthal and polar angles



% set initial condition
positions(1, :) = initial_position;
velocities(1, :) = initial_velocity;

% calculate Lorentz force with numerical integration!
for step = 1:num_steps
    % Calculate Lorentz force
    F = q * cross(velocities(step, :), B);
    
    % Update velocity and position using the equations of motion
    acceleration = F / m;
    velocities(step + 1, :) = velocities(step, :) + acceleration * dt;
    positions(step + 1, :) = positions(step, :) + velocities(step + 1, :) * dt;

    % Convert velocity vector to spherical coordinates
    [azimuthal_angle, polar_angle, ~] = cart2sph(velocities(step + 1, 1), velocities(step + 1, 2), velocities(step + 1, 3));

    angles(step + 1, :) = [azimuthal_angle, polar_angle];
end

% Calculate initial and final angles
initial_angles = angles(1, :);
final_angles = angles(end, :);

% Display results
figure;
plot3(positions(:, 1), positions(:, 2), positions(:, 3));
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
title('Proton Trajectory in Magnetic Field');
grid on;

% Display initial and final angles
disp('Initial Angles (Azimuthal, Polar):');
disp(rad2deg(initial_angles));
disp('Final Angles (Azimuthal, Polar):');
disp(rad2deg(final_angles));

% Display end position
disp('End Position:');
disp(positions(end, :));