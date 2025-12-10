% Constants (needed for velocity calculation for example)
atomic_mass_unit_MeV_c2 = 931.494; % MeV/c^2
c = 299792458 * 100; % [mm/s]

% Parameters SPR_values.csv watertestSPR.csv
SPR_values = readmatrix('watertestSPR.csv');
E0_MeV = 100; % Initial proton energy
B = [0; 0; 3]; % Magnetic field [T]
grid_step = 0.00109375; % [m] This is the CT grid step /"PixelSpacing"  one can aquire by looking at metadata using 3D Slicer for example
%total_distance = 83.91796875; % [mm] CT rows * CT columns * CT grid step
total_distance = 2000;

% Possible needed pre-calculations
gamma0 = 1 + E0_MeV / atomic_mass_unit_MeV_c2;
beta0 = sqrt(1 - 1/gamma0^2);
v0 = beta0 * c;

% Parameters - Initialization
initial_position = [0; 0; 0];
%initial_position = [0; 11.2; 0]; % Initial position [mm] Current analyzed CT image has 225 rows -> 225/2=112.5. For our simulation we need 11.25 as an input for the beam to start at the center.
initial_velocity = [v0; 0; 0]; % Initial velocity [mm/s] Here in x-direction

% Create an instance of ProtonSimulation
vacuumSim = Vacuum_Simulation(SPR_values, E0_MeV, B, grid_step, total_distance, initial_position, initial_velocity);

% Initialize and run the simulation
vacuumSim = vacuumSim.initializeStep();
vacuumSim = vacuumSim.simulate();
vacuumSim.saveResults();

% Display the number of steps with non-zero values
vacuumSim.displayStepRange(); % You can input step values up until this output value!

% Accessing specific values
step = 7200;
disp(['Energy at step ', num2str(step), ': ', num2str(vacuumSim.getEnergyAtStep(step)), ' MeV']);
disp(['X-position at step ', num2str(step), ': ', num2str(vacuumSim.positions(1, step)), ' cm']);
