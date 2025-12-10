% RAMDIM: Raytracing Algorithm for Magnetic Deflection of Ions in Media

% Define constants
E0 = 150; % Initial energy of the proton beam in MeV
B0 = 3; % Magnetic flux density in Tesla
d_air = 25; % Thickness of the air gap in cm
alpha = 2.43e-3; % Constant alpha in MeV^(-p)cm
p = 1.75; % Power-law exponent
q = 1.602e-19; % Elementary charge in Coulombs
m = 1.67262192e-27; % Mass of a proton in kg
c = 299792458; % Speed of light in m/s

% Calculate the range in water
R0 = alpha * E0^p;

% Calculate energy from MeV to J for further calculations
E_SI = E0*10^6 * 1.602176634 * 10^(-19);
disp(['with given energy ',num2str(E0), ' MeV, we get Joule ', num2str(E_SI)]);




% Calculate the non-relativistic initial velocity
v_initial = sqrt(2*E_SI/m);
disp(['the inital non relativistic velocity is ', num2str(v_initial)]);

% Calculate the lorentz factor and the relativistic velocity
gamma_initial = 1 / sqrt(1 - (v_initial^2 / c^2)); % Lorentz factor
disp(['the relativistic lorentz factor is ', num2str(gamma_initial)]);

v_rel = gamma_initial * v_initial; % Relativistic initial velocity of the proton beam
disp(['the relativistic velocity is ', num2str(v_rel)]);

% Calculate initial gyroradius by magnetic field
r_initial = (m * v_initial * gamma_initial) / (q * B0); 



% Define the step size for the trajectory calculation
ds = 0.1; % Step size in cm

% Define the number of steps for the trajectory calculation
n = round(R0 / ds);

% Initialize the proton position and velocity
x = [0; 0; 0]; % Position in cm
v = [v_rel; 0; 0]; % Velocity in m/s

%% Initialize the traveled path length in water
%s_water = 0;

% Initialize arrays to store the positions
x_positions = zeros(3, n+1); % Array to store x positions
s_water = zeros(1, n+1); % Array to store the traveled path length in water

% Store the initial position
x_positions(:,1) = x;
s_water(1) = 0;

% Perform the trajectory calculation using the RAMDIM model
for i = 1:n
    % Calculate the traveled path length in water at each step
    Ei = E0 - (q * B0 * ds / (2 * m * c)) * (2 * i - 1); % Calculate the energy at step i
    si = R0 - alpha * Ei^p;
    
    % Accumulate the traveled path length in water
    s_water = s_water + si;
    
    % Calculate the gyroradius and deflection angle
    r = m * v(1) / (q * B0); %%%%% Velocity should be calculated new every time??!!
    phi = ds / r;
    
    % Calculate the new position and velocity
    x = x + ds * v / 100; % Convert velocity to cmDCDw
    v = [v(1) * cos(phi); v(1) * sin(phi); v(3)];

    % Store the new position
    x_positions(:,i+1) = x;
    s_water(i+1) = s_water(i) + si;
end

% Display the total traveled path length in water
disp(['Total traveled path length in water: ', num2str(s_water), ' cm']);

% Plot the trajectory
figure;
plot3(x_positions(1,:), x_positions(2,:), x_positions(3,:), '-o');
xlabel('X position (cm)');
ylabel('Y position (cm)');
zlabel('Z position (cm)');
title('Proton Beam Trajectory');
grid on;