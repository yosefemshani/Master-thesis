% Constants
m = 1.6726219e-27; % Mass of proton in kg
q = 1.60217662e-19; % Charge of proton in Coulombs
c = 2.998e8; % Speed of light in m/s
v0 = c * sqrt(1 - (m * c^2 / (200 * 10^6 * q))^2); % Initial velocity in m/s

% Magnetic field components
B = [0; 0; 0.5]; % Magnetic field vector in Tesla

% Initial conditions
r0 = [0; 0; 0]; % Initial position vector in meters
v0_vector = [0; 0; v0]; % Initial velocity vector in m/s

% Time parameters
t_start = 0;
t_end = 10; % You can change this depending on the desired time range
dt = 1e-10; % Time step size

% Number of time steps
num_steps = ceil((t_end - t_start) / dt);

% Preallocate array for storing trajectory
r = zeros(3, num_steps); % Matrix to store position vectors

% Initial values
r(:,1) = r0;

% Main loop to calculate trajectory
for i = 2:num_steps
    % Calculate the Lorentz force components
    F = q * cross(v0_vector, B);
    
    % Update accelerations using F = ma (Lorentz force)
    a = F / m;
    
    % Update positions using v = v0 + at
    r(:,i) = r(:,i-1) + v0_vector * dt;
    
    % Update velocities using a = F/m
    v0_vector = v0_vector + a * dt;
end

% Plot the trajectory
figure;
plot3(r(1,:), r(2,:), r(3,:));
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('Trajectory of Proton under Magnetic Field');
grid on;
