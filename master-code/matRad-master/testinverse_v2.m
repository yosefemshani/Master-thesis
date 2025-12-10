% Inverse dynamics simulation

% definition of basic 3D 10x10x10 cm^3 box
x = 10;
y = 10;
z = 10;

% definition of parameters
q = 1.602 * 10^(-19); % [C]
m = 1.6726 * 10^(-27); % [kg]
B = [0, 0, 1.5]; % 1.5 [T] in z-axis
E0 = 150; % MeV
c = 299792458; % Speed of light in m/s
d_air = 25; % Define Thickness of the air gap [cm] (POSSIBLE FURTHER CONSIDERATIONS NEEDED!?)


% Calculate energy from MeV to J for further calculations
E_SI = E0*10^6 * 1.602176634 * 10^(-19);

% Calculate the initial lorentz factor and the relativistic velocity []
gamma_initial = 1 / sqrt(1 - (v_initial^2 / c^2)); % Lorentz factor

%initial velocity
v_initial = sqrt(2*E_SI/m);

% Calculate the initial relativistic velocity [m/s]
v_rel = gamma_initial * v_initial;

%initial gyroradius
r_initial = (m * v_initial * gamma_initial) / (q * B0(3)); 

% Rotation angle (in radians) for counter-clockwise rotation around the z-axis
rotation_angle = asin(d_air/r_initial); % Adjust the angle as needed

% Rotation matrix for counter-clockwise rotation around the z-axis
R = [cos(rotation_angle), -sin(rotation_angle), 0;
     sin(rotation_angle), cos(rotation_angle), 0;
     0, 0, 1];



% initial conditions
initial_position = [0,0,0]; % initial position of x y and z in [m]
rotated_initial_position = R * initial_position'; 
%initial_velocity = [269813212.2, 0, 0]; % initial velocity v = 0.9c in x-axis [m/s]
initial_velocity = [v_rel, 0, 0]; % initial velocity v = 0.9c in x-axis [m/s]


% Desired end position and angles
desired_end_position = [1.0e+04 * (-5.1529), 0, 1.0e+04 * 0.5570];
desired_end_angles = deg2rad([180, -75.6570]); % Set your desired final angles

% Simulation parameters
tolerance = 1e-5; % Tolerance for position convergence
max_iterations = 1000; % Maximum number of iterations

% Initialize arrays to store position, velocity, and angles at each step
positions = zeros(max_iterations + 1, 3);
velocities = zeros(max_iterations + 1, 3);
angles = zeros(max_iterations + 1, 2);

% Set initial conditions
positions(1, :) = rotated_initial_position';
velocities(1, :) = initial_velocity;

% Perform inverse dynamics simulation
for iteration = 1:max_iterations
    % Calculate Lorentz force
    F = q * cross(velocities(iteration, :), B);
    
    % Update velocity and position using the equations of motion
    acceleration = F / m;
    velocities(iteration + 1, :) = velocities(iteration, :) + acceleration * dt;
    positions(iteration + 1, :) = positions(iteration, :) + velocities(iteration + 1, :) * dt;

    % Convert current position to spherical coordinates
    [azimuthal_angle, polar_angle, ~] = cart2sph(positions(iteration + 1, 1), positions(iteration + 1, 2), positions(iteration + 1, 3));

    % Store current angles
    angles(iteration + 1, :) = [azimuthal_angle, polar_angle];
    
    % Check for convergence
    if norm(positions(iteration + 1, :) - desired_end_position) < tolerance && all(abs(angles(iteration + 1, :) - desired_end_angles) < tolerance)
        disp(['Converged after ', num2str(iteration), ' iterations.']);
        break;
    end
end

% Display results
figure;
plot3(positions(:, 1), positions(:, 2), positions(:, 3));
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
title('Inverse Dynamics: Proton Trajectory in Magnetic Field');
grid on;

% Display initial and final angles
disp('Desired End Angles (Azimuthal, Polar):');
disp(rad2deg(desired_end_angles));
disp('Final Angles (Azimuthal, Polar):');
disp(rad2deg(angles(end, :)));

% Display end position
disp('Desired End Position:');
disp(desired_end_position);
disp('Final Position:');
disp(positions(end, :));