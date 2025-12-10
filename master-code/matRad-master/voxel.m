%% Integrating magnetic field grid size model with consideration of proton energy loss applied on a given CT patient data set (2_04_P_real_cropped)
%% Author: Yosef Emshani, 01.03.2024 17:51


%% Clear workspace

clear;  % Clear workspace
clc;    % Clear command window
close all;  % Close all figures


%% Definition of constant variables

mp_MeV = 938.272013;
mp = 1.67262158e-27; % [kg]
c = 299792458*100; % [mm/s]
atomic_mass_unit_MeV_c2 = 931.494; %MeV/c^2; double check this value later once code runs! might be 931.494e6... 
q = 1.602e-19; % elementary charge: positive value for protons!


%% Definition of variables one might change

E0_MeV = 100; % initial proton energy, 332 MeV drops to 0, 333 doesnt..
B = [0; 0; 1.5]; % [T]


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

%% Import the SPR values from the CSV file
SPR_values = readmatrix('IMG0035_cropped_SPR.csv');
[num_rows, num_cols] = size(SPR_values);

% Ensure SPR_values dimensions match CT data
assert(numel(SPR_values) == num_rows * num_cols, 'SPR values dimensions do not match CT data dimensions');

% Initialize array to store energy losses adjusted by SPR
energy_losses_adjusted = zeros(1, num_steps);

% Loop through each step
for step = 2:num_steps
    % Calculate time taken to travel each step
    time_step = grid_step / norm(velocities(:, step - 1));

    F = q * cross((velocities(:, step - 1)), B);
    a = F/mp;
    
    % Calculate energy loss without SPR adjustment
     energy_loss_keV_um = calculate_energy_loss(E_MeV);

    % Check if energy is below the BETHE_LOWER_LIMIT
    if E_MeV < 0.49
        % Set energy loss to 0
        energy_loss_keV_um = 0;
        % Set energy to 0
        E_MeV = 0;
        % Break out of the loop
        disp(['Energy dropped below Bethe lower limit, E=0 MeV']);
        break;
    else
        % Find the corresponding SPR value using interpolation
        x_pos = positions(1, step-1) * 10; % Convert to correct scale
        y_pos = positions(2, step-1) * 10; % Convert to correct scale

        % Ensure positions are within the CT grid boundaries
        x_index = min(max(round(x_pos), 1), num_cols);
        y_index = min(max(round(y_pos), 1), num_rows);

        % Debug
        disp(['x index is at ', num2str(x_index), ' [idk]']);
        disp(['y index is at ', num2str(y_index), ' [idk]']);

        % Get SPR value for the current position
        SPR_value = SPR_values(y_index, x_index);

        % Adjust energy loss by multiplying with SPR value for this step
        energy_loss_keV_um_adjusted =  energy_loss_keV_um * SPR_value;

        % Store adjusted energy loss
        energy_losses_adjusted(step) =  energy_loss_keV_um_adjusted;

        % Update energy for the next step
        E_MeV = E_MeV -   energy_loss_keV_um_adjusted * grid_step;
    end

    % Recalculate gamma and beta
    gamma = 1 + E_MeV / atomic_mass_unit_MeV_c2;
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
    y_new = positions(2, step - 1) + vy_new * (grid_step / abs(vx_new));
    z_new = positions(3, step - 1) + vz_new * (grid_step / abs(vx_new));
    
    % Update positions
    positions(:, step) = [x_new; y_new; z_new];
    
    % Update remaining distance after the step
    x_distance = 0;
end

    % Store time needed for x-y position to get to grid step
    times(step) = times(step - 1) + time_step;
    
    % Store accelerations
    accelerations(:, step) = a;

end

% Create a logical mask to identify rows that do not contain all zeros
non_zero_mask = any(positions ~= 0, 1);

% Filter the positions array to only include non-zero rows
filtered_positions = positions(:, non_zero_mask);

% Save filtered positions to CSV file
csvwrite('trajectory.csv', filtered_positions');
disp(['Non-zero positions have been saved!']);

% Plot trajectory
figure;
plot(filtered_positions(1, :), filtered_positions(2, :), 'b.-', 'LineWidth', 1);
title('Trajectory of Proton under Magnetic Field with Energy Loss');
xlabel('X Position (cm)');
ylabel('Y Position (cm)');
grid on;
axis equal;


function energy_loss_MeV_cm = calculate_energy_loss(E_MeV)
    % constant variables
    me_MeV = 0.510998918; %(see https://github.com/libamtrack/library/blob/master/include/AT_Constants.h#L38 l. 38)
    mp_MeV = 938.272013; %(see https://github.com/libamtrack/library/blob/master/include/AT_Constants.h#L35 l. 35)
    m_MeV = 1.0079 * mp_MeV; %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L278 l. 280)
    c = 299792458; % Speed of light in m/s
    re = 2.8179403262 * 10^(-15); % electron radius in m
    NA = 6.02214086 * 10^23; % 1/mol
    atomic_mass_unit_MeV_c2 = 931.494; %MeV/c^2
    Dirac_constant_J_s = 1.054571628e-34;
    BETHE_LOWER_LIMIT_E_MEV_U = 0.49;
    phase_undefined = 0;
    phase_condensed = 1;
    phase_gas = 2;
    Z = 1; %for protons

    % variables one might change!
    E_restricted_keV = 1000;
    I_eV = 78; % 75 eV for liquid water. for others, see https://github.com/libamtrack/library/blob/44dc48cfa977c05008ad12646d798d6c4b6ea504/include/AT_DataMaterial.h l. 197-203
    AT = 18.0006; % mass number for H2O
    ZT = 10; % ordinary number for target: water (H2O: 2 * 1 + 1 * 8 = 10)
    phase = phase_condensed; %for liquid water! for others, see https://github.com/libamtrack/library/blob/master/include/AT_DataMaterial.h#L246
                                                               %lines 240-253 for different materials!

    % simple calculations and conversion factors
    gamma = 1 + E_MeV / atomic_mass_unit_MeV_c2; %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L45 lines 65-82)
    assert(gamma >= 1.0);

    beta = sqrt(1 - 1/gamma^2); %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L45 line 75)
    mass_correction_term = 1 + (2 * (me_MeV / m_MeV) / gamma) + (me_MeV / m_MeV)^2; %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L278 l. 281)
    Wm_MeV = 2 * me_MeV * beta^2 / (1 - beta^2) / mass_correction_term; %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L285 l.288)
    beta2 = beta^2;
    I_MeV = I_eV * 1e-6;
    m_to_cm = 100;
    MeV_to_J = 1.60217646e-13;
    I_J = I_MeV * MeV_to_J;

    % Stopping power

    % Restricted stopping number requested?
    restricted = false;
    if E_restricted_keV > 0.0 && (E_restricted_keV / 1000.0) < Wm_MeV
        restricted = true;
    end

    % First part of stopping number
    SN11 = (2.0 * me_MeV * beta2) / (1.0 - beta2);
    assert(I_MeV > 0);
    SN11 = SN11 / I_MeV;

    if restricted
        Wm_MeV = E_restricted_keV * 1e-3;
    end
    SN12 = Wm_MeV / I_MeV;

    % Second part of stopping number
    SN2 = beta2;
    if restricted
        SN2 = SN2 / 2;
        SN2 = (SN2 + (1.0 - beta2) * Wm_MeV) / (4.0 * me_MeV);
    end

    % Third part of stopping number (density correction following Sternheimer, 1971)
    delta = 0.0;
    if phase ~= phase_undefined

        %MIGHT NEED TO VERIFY IF THIS KINETIC_VARIBALE & PLASMA ENERGY CALC IS CORRECT!
        gamma_single =  1.0 + E_MeV / atomic_mass_unit_MeV_c2;
        assert(gamma_single >= 1.0);

        assert(E_MeV >= 0);
        beta_single = sqrt(1 - (1/gamma_single^2));

        assert(beta_single * gamma_single > 0);
        kinetic_variable_single = log10(beta_single*gamma_single);
        kinetic_variable = kinetic_variable_single;

        rho_gcm3 = 1;
        numbers_of_atoms_per_g = NA / AT;
        numbers_of_electrons_per_g = numbers_of_atoms_per_g * ZT;
        electron_density_per_cm3 = numbers_of_electrons_per_g * rho_gcm3;
        electron_density_m3_single = electron_density_per_cm3 * m_to_cm * m_to_cm * m_to_cm;
        electron_density_m3 = electron_density_m3_single;

        %https://github.com/libamtrack/library/blob/master/src/AT_DataMaterial.c#L294
        %line 307
        plasma_energy_J = sqrt((4 * pi * electron_density_m3 * re) * Dirac_constant_J_s * c);

        C = 1.0 + 2.0 * log(I_J / plasma_energy_J);

        % Find x_0 and x_1 dependent on phase, I-value, and C
        if phase == phase_condensed
            if I_eV < 100
                x_1 = 2.0;
                if C <= 3.681
                    x_0 = 0.2;
                else
                    x_0 = 0.326 * C - 1.0;
                end
            else % I_eV >= 100
                x_1 = 3.0;
                if C <= 5.215
                    x_0 = 0.2;
                else
                    x_0 = 0.326 * C - 1.5;
                end
            end
        else % gaseous material
            x_0 = 0.326 * C - 2.5;
            x_1 = 5.0;
            if C < 10.0
                x_0 = 1.6;
                x_1 = 4.0;
            end
            if C >= 10.0 && C < 10.5
                x_0 = 1.7;
                x_1 = 4.0;
            end
            if C >= 10.5 && C < 11.0
                x_0 = 1.8;
                x_1 = 4.0;
            end
            if C >= 11.0 && C < 11.5
                x_0 = 1.9;
                x_1 = 4.0;
            end
            if C >= 11.5 && C < 12.25
                x_0 = 2.0;
                x_1 = 4.0;
            end
            if C >= 12.25 && C < 13.804
                x_0 = 2.0;
                x_1 = 5.0;
            end
        end

        x_a = C / 4.606;
        m = 3.0;
        a = 4.606 * (x_a - x_0) / ((x_1 - x_0)^m);

        if kinetic_variable >= x_0 && kinetic_variable <= x_1
            delta = 4.606 * kinetic_variable - C + a * (x_1 - kinetic_variable)^m;
        end
        if kinetic_variable > x_1
            delta = 4.606 * kinetic_variable - C;
        end
    end
    SN3 = delta;

    % Forth part of stopping number (shell correction) TODO: implement

    %assert(SN11 > 0);
    %assert(SN12 > 0);

    stopping_number = 0.5 * log(SN11 * SN12) - SN2 - SN3;

    % Leading energy loss term
    %(we assume that effective charge is used! otherwise, change code with the help of https://github.com/libamtrack/library/blob/44dc48cfa977c05008ad12646d798d6c4b6ea504/src/AT_StoppingPowerDataBethe.c#L49 l. 62-65)

    assert(AT > 0);
    assert(beta2 > 0);

    % calculation of effective charge according to Barkas-Bethe approximation
    if Z ~= 1
        effective_charge = Z * (1.0 - exp(-125.0 * beta / (Z^(2.0/3.0))));
    else
        effective_charge = 1.0 - exp(-125.0 * beta);
    end

    %since we only use effective charge we can say z = Z;
    z = Z;

    % ICRU49, p.6, after Cohen and Taylor (1986), k_MeV_cm2_g = 0.307075
    energy_loss_leading_term_MeV_cm2_g = 0.307075 * (ZT / AT) * (z^2) / beta2;

    
    % Energy loss
    % Compute only above 1.0 MeV, otherwise theory is too wrong below return zero
    % TODO: Find smarter criterion because this may cause problems in the code (as it did
    % TODO: with the inappropriately set lower limit for CSDA range integration (was 0, now 1.0 MeV)
    
    if E_MeV >= BETHE_LOWER_LIMIT_E_MEV_U
        energy_loss_MeV_cm = energy_loss_leading_term_MeV_cm2_g * stopping_number;
    else
       energy_loss_MeV_cm = 0;
    end

    % unit conversion
    energy_loss_keV_um = energy_loss_MeV_cm / 10; % [keV/μm] AND [MeV/mm]
    energy_loss_MeV_m = energy_loss_MeV_cm / 100; % [MeV/m]

    % TEST OUTPUT DISPLAY FOR CANDIDATE VARIABLES
    disp(['dE/dx is ', num2str(energy_loss_keV_um), ' [MeV/mm]']);
    disp(['current Energy is ', num2str(E_MeV), ' [MeV]']);
    disp(['---------------------------------------------------']);

    %disp(['dE/dx is ', num2str(energy_loss_keV_um), ' [keV/μm]']);
end