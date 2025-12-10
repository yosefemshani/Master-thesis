% run_ProtonSimulation.m

% Accessing specific values
%step = 7200;
%disp(['Energy at step ', num2str(step), ': ', num2str(protonSim.getEnergyAtStep(step)), ' MeV']);
%disp(['X-position at step ', num2str(step), ': ', num2str(protonSim.positions(1, step)), ' mm']);




% Example usage of the new function
y_position = 11.2; % Example initial y-position
%y_position = 14.239888751946; % Example initial y-position
E0_MeV = 225; %159.859148730345; % Example initial energy in MeV
magnetic_field = [0; 0; 0]; % Example magnetic field in T
SPR = readmatrix('watertestSPR.csv'); 
% Load SPR values SPR_values(_slice).csv watertestSPR.csv bone_SPR.csv heavybone_SPR.csv




final_pos = analyzeProtonTrajectory(y_position, E0_MeV, magnetic_field, SPR);
disp([num2str(final_pos(1)')]);




%%%%%%%%%%%%%%%%%% Liheng write here
%
%
%Target = analyzeProtonTrajectory(y_position, E0_MeV, [0;0;0]);
%Target_start = analyzeProtonTrajectory(y_position, E0_MeV, magnetic_field);
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
%    position_this_step = analyzeProtonTrajectory(yi, Ei, magnetic_field);
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





