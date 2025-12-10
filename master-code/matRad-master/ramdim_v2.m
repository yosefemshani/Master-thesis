% RAMDIM: Raytracing Algorithm for Magnetic Deflection of Ions in Media
% source: https://doi.org/10.1088/1361-6560/62/4/1548
% Author of this algorithm implementation: Yosef Emshani, 23.1.2024 15:08

% Define constants
E0 = 150; % Initial energy of the proton beam in MeV
B0 = [0, 0, 1.5]; % Magnetic flux density in Tesla
alpha = 2.43e-3; % Constant alpha in MeV^(-p)cm
p = 1.75; % Power-law exponent
q = 1.602e-19; % Elementary charge in Coulombs
m = 1.67262192e-27; % Mass of a proton in kg
c = 299792458; % Speed of light in m/s


% Calculate energy from MeV to J for further calculations
E_SI = E0*10^6 * 1.602176634 * 10^(-19);
disp(['with given energy ',num2str(E0), ' MeV, we get Joule ', num2str(E_SI)]);



% Calculate the non-relativistic initial velocity [m/s]
v_initial = sqrt(2*E_SI/m);
disp(['the inital non relativistic velocity is ', num2str(v_initial), ' m/s']);

% Calculate the initial lorentz factor and the relativistic velocity []
gamma_initial = 1 / sqrt(1 - (v_initial^2 / c^2)); % Lorentz factor
disp(['the relativistic lorentz factor is ', num2str(gamma_initial)]);


% Calculate the initial relativistic velocity [m/s]
v_rel = gamma_initial * v_initial;
disp(['the relativistic velocity is ', num2str(v_rel), ' m/s']);

% Calculate initial gyroradius by magnetic field  [m]
r_initial = (m * v_initial * gamma_initial) / (q * B0(3)); 
disp(['the initial gyroradius is  ', num2str(r_initial), ' m']);



% Calculate the range in water [cm]
R0 = alpha * E0^p;
disp(['the range in water is  ', num2str(R0), ' cm']);

% Define Thickness of the air gap [cm] (POSSIBLE FURTHER CONSIDERATIONS NEEDED!?)
d_air = 25;


% Define the step size for the trajectory calculation
ds = 0.1; % Step size in cm

% Define the number of steps for the trajectory calculation
n = round(R0 / ds);

% Initialize the proton position and velocity
x = [0; 0; 0]; % Initial position in cm
v = [v_rel; 0; 0]; % Initial velocity in m/s


%% Understand that there is a difference between the entrance position
%% of the proton beam at the surface of the water phantom and
%% the following positions inside the water phantom! Ref. eq. 16


% Define initial angle
angle_initial = asin(d_air/r_initial);

% Define 3x3 rotation matrix which rotates counter clockwise around z-axis
Rotation = [cos(angle_initial) -sin(angle_initial) 0;
            sin(angle_initial)  cos(angle_initial) 0;
            0           0          1];

% Calculate entrance position x1 of the proton beam at the surface of the
% water phantom
x1 = Rotation * x;
disp(['the entrance position on the surface is  [', num2str(x1(1)), ', ', num2str(x1(2)), ', ', num2str(x1(3)), ']']);

%Define time steps [s]
dt = 1e-9;

% Calculate the velocity coordinates
vx = v_initial * (cos((q * B0(3)*dt) / (gamma_initial * m)));
vy = v_initial * (sin((q * B0(3)*dt) / (gamma_initial * m)));
vz = 0;



%% Tracjectory calculation

% Initialize arrays to store position, velocity, and angles at each step
positions = zeros(n + 1, 3);
velocities = zeros(n + 1, 3);

% Setup array for initial conditions
%positions = [0,0,0];
%velocities = [v_rel,0,0];

positions(1, :) = x1;
velocities(1, :) = [v_rel,0,0];



% calculate Lorentz force with numerical integration!
for step = 1:n
    % Calculate Lorentz force
    F = q * cross(velocities(step, :), B0);
    
    % Update velocity and position using the equations of motion
    acceleration = F / m;

    velocities(step + 1, :) = velocities(step, :) + acceleration * dt;
    positions(step + 1, :) = positions(step, :) + velocities(step + 1, :) * dt;
end

disp([num2str(positions)]);

%function F = calculateLorentzForce(v, B0, q)
%    % Calculates the Lorentz force based on velocity v, magnetic field B0, and charge q
%    F = q * cross(v, B0);
%end


% Display results
figure;
plot3(positions(:, 1), positions(:, 2), positions(:, 3));
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
zlabel('Z-axis (m)');
title('Proton Trajectory in Magnetic Field');
grid on;

%% Energy loss calculation using Bethe bloch formula


% Define material properties for liquid water
% Values taken from NIST: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
Z = 10;         % Atomic number for water (approximately)
A = 18.01528;   % Atomic mass for water (approximately)
I = 75.0;       % Mean excitation potential for water (eV)
ZoverA = Z/A;   % Z/A ratio

% Define particle properties
E = 1e6;        % Particle energy in MeV (for example, 1 MeV)
z = 1;          % Particle charge (for example, for a proton)

% Calculate stopping power
[dEdx, ShellCorr, DensityCorr] = BetheBloch(E, Z, A, ZoverA, I, z);

% Display results
fprintf('Stopping Power: %f MeV/cm\n', dEdx);
fprintf('Shell Correction: %f MeV/cm\n', ShellCorr);
fprintf('Density Correction: %f MeV/cm\n', DensityCorr);


function [dEdx, ShellCorr, DensityCorr] = BetheBloch(E, Z, A, ZoverA, I, z)
    % Constants
    me = 0.511; % MeV/c^2, electron mass
    K = 0.307; % MeV cm^2/mol, constant
    M = 9.10938356e-31; % kg, electron mass
    e = 1.60217662e-19; % C, elementary charge
    Na = 6.022e23; % mol^-1, Avogadro's number

    % Parameters
    beta = sqrt(E.^2./(E.^2 + 2*E*me)); % Velocity
    gamma = 1./sqrt(1 - beta.^2); % Lorentz factor
    Tmax = 2*me.*beta.^2.*gamma.^2./(1 + 2*gamma.*me./M + (me./M).^2); % Maximum kinetic energy transfer
    
    % Bethe-Bloch formula
    dEdx = K * ZoverA ./ beta.^2 .* (0.5 .* log(2*me.*beta.^2.*gamma.^2.*Tmax./(I^2)) - beta.^2 - delta(Z, beta, me, I)); 
    
    % Shell correction
    ShellCorr = log(2*me*beta.^2*gamma.^2.*Tmax./(I^2)) - 2 .* beta.^2;
    
    % Density correction
    X0 = 36.08; % cm, radiation length for air
    m = 3; % exponent
    C = 5.2146; % empirical constant for density correction
    X = log10(beta.*gamma);
    a = log10(X0);
    DensityCorr = log((X.^m + C)./(a.^m + C));
    
    dEdx = dEdx - ShellCorr - DensityCorr;
end


function d = delta(Z, beta, me, I)
    a = 0.5 * log(2 * me * beta.^2);
    b = log(I);
    z = log(Z);
    d = a - b - z;
end

