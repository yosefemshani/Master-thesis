% Constants (needed for velocity calculation for example)
atomic_mass_unit_MeV_c2 = 931.494; % MeV/c^2
c = 299792458 * 100; % [mm/s]

% Parameters
SPR_values = readmatrix('watertestSPR.csv');
initial_position = [0; 11.2; 0]; % Starting position [mm]
end_position = [0; 20; 0]; % Desired end-position [mm]
B = [0; 0; 1.5]; % Magnetic field [T]
grid_step = 0.00109375; % [m]
total_distance = 83.91796875; % [mm]
initial_velocity = [0; 0; 0]; % Assuming initial velocity is zero for simplicity

% Create an instance of ProtonSimulation_inverse
protonSim = ProtonSimulation_inverse(SPR_values, end_position, B, grid_step, total_distance, initial_position, initial_velocity);

% Find the required energy to reach the end position
required_energy_MeV = protonSim.findRequiredEnergy();

% Display the required energy
disp(['Required energy to reach the end position: ', num2str(required_energy_MeV), ' MeV']);

% Save the results
protonSim.saveResults();
