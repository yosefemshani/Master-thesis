classdef HeliumSimulation
    properties
        % Constant variables
        mp_MeV = 3727.379; % Mass of helium nucleus [MeV/c²]
        mp = 6.644657230e-27; % Mass of helium nucleus [kg]
        c = 299792458 * 100; % [mm/s]
        atomic_mass_unit_MeV_c2 = 931.494; % MeV/c²
        q = 2 * 1.602e-19; % Elementary charge for helium (2 protons)

        % Variables one might change
        E0_MeV; % Initial helium energy
        B; % Magnetic field [T]
        initial_position; % Initial position [mm]
        initial_velocity; % Initial velocity [mm/s]
        
        % Arrays to store results
        positions;
        velocities;
        accelerations;
        times;
        energy_losses_adjusted;
        energies; % To store energy at each step
        
        % Initial conditions
        grid_step; % [m]
        total_distance; % [mm]
        num_steps;
        SPR_values;
        E_MeV;
        x_distance = 0;
    end
    
    methods
        function obj = HeliumSimulation(SPR_values, E0_MeV, B, grid_step, total_distance, initial_position, initial_velocity)
            % Constructor to initialize properties and load SPR values
            obj.SPR_values = SPR_values;
            obj.E0_MeV = E0_MeV;
            obj.B = B;
            obj.grid_step = grid_step;
            obj.total_distance = total_distance;
            obj.initial_position = initial_position;
            obj.initial_velocity = initial_velocity;
            steps = round(obj.total_distance / obj.grid_step);
            obj.num_steps = steps;
            obj.positions = zeros(3, steps);
            obj.velocities = zeros(3, steps);
            obj.accelerations = zeros(3, steps);
            obj.times = zeros(1, steps);
            obj.energy_losses_adjusted = zeros(1, steps);
            obj.energies = zeros(1, steps); % Initialize energies array
            obj.E_MeV = obj.E0_MeV;
        end
        
        function obj = initializeStep(obj)
            % Initialize first step of velocity and acceleration
            gamma0 = 1 + obj.E0_MeV / obj.atomic_mass_unit_MeV_c2;
            beta0 = sqrt(1 - 1/gamma0^2);
            v0 = beta0 * obj.c;
            disp(['Initial velocity is : ', num2str(v0)]);
            v = obj.initial_velocity; % Use the initial velocity input
            F = - obj.q * cross(v, obj.B);
            a = F / obj.mp;
            obj.positions(:, 1) = obj.initial_position; % Use the initial position input
            obj.velocities(:, 1) = v;
            obj.energies(1) = obj.E0_MeV; % Store initial energy
        end
        
        function obj = simulate(obj)
            [num_rows, num_cols] = size(obj.SPR_values);
            
            for step = 2:obj.num_steps
                % Time step calculation based on half the grid step
                half_grid_step = obj.grid_step / 2;
                half_time_step = half_grid_step / norm(obj.velocities(:, step - 1));

                % First half-step to mid-point of the voxel
                F_mid = - obj.q * cross(obj.velocities(:, step - 1), obj.B);
                a_mid = F_mid / obj.mp;

                mid_velocities = obj.velocities(:, step - 1) + a_mid * half_time_step;
                mid_positions = obj.positions(:, step - 1) + obj.velocities(:, step - 1) * half_time_step;

                % Calculate energy loss at the mid-point
                energy_loss_keV_um = obj.calculate_energy_loss(obj.E_MeV);

                if obj.E_MeV < 0.49
                    energy_loss_keV_um = 0;
                    obj.E_MeV = 0;
                    disp(['Energy dropped below Bethe lower limit, E=0 MeV']);
                    break;
                else
                    x_pos = mid_positions(1) * 10; % Convert to mm
                    y_pos = mid_positions(2) * 10; % Convert to mm
                    x_index = min(max(round(x_pos), 1), num_cols);
                    y_index = min(max(round(y_pos), 1), num_rows);
                    SPR_value = obj.SPR_values(y_index, x_index);
                    energy_loss_keV_um_adjusted = energy_loss_keV_um * SPR_value;
                    obj.energy_losses_adjusted(step) = energy_loss_keV_um_adjusted;
                    obj.E_MeV = obj.E_MeV - energy_loss_keV_um_adjusted * obj.grid_step;
                    obj.energies(step) = obj.E_MeV; % Store energy at current step
                end

                % Update gamma and beta
                gamma = 1 + obj.E_MeV / obj.atomic_mass_unit_MeV_c2;
                beta = sqrt(1 - 1/gamma^2);
                if gamma < 1
                    disp(['Gamma dropped below 1 at step ', num2str(step)]);
                    break;
                end

                % Complete the second half-step to the end of the voxel
                F_end = - obj.q * cross(mid_velocities, obj.B);
                a_end = F_end / obj.mp;

                obj.velocities(:, step) = mid_velocities + a_end * half_time_step;
                obj.positions(:, step) = mid_positions + mid_velocities * half_time_step;
                obj.accelerations(:, step) = a_end;
                obj.times(step) = obj.times(step - 1) + 2 * half_time_step;
            end
        end
        
        function energy_loss_MeV_cm = calculate_energy_loss(obj, E_MeV)
            % constant variables
            me_MeV = 0.510998918; % Electron mass [MeV/c²]
            mp_MeV = 3727.379; % Helium mass [MeV/c²]
            m_MeV = 4 * mp_MeV; % Mass of helium nucleus
            c = 299792458; % Speed of light in m/s
            re = 2.8179403262 * 10^(-15); % electron radius in m
            NA = 6.02214086 * 10^23; % 1/mol
            atomic_mass_unit_MeV_c2 = 931.494; % MeV/c²
            Dirac_constant_J_s = 1.054571628e-34;
            BETHE_LOWER_LIMIT_E_MEV_U = 0.49;
            phase_undefined = 0;
            phase_condensed = 1;
            phase_gas = 2;
            Z = 2; % Atomic number for helium

            % variables one might change!
            E_restricted_keV = 1000;
            I_eV = 78; % Ionization energy for water
            AT = 18.0006; % Mass number for H2O
            ZT = 10; % Atomic number for target: water (H2O)
            phase = phase_condensed; % Phase of the material

            % simple calculations and conversion factors
            gamma = 1 + E_MeV / atomic_mass_unit_MeV_c2;
            assert(gamma >= 1.0);

            beta = sqrt(1 - 1/gamma^2);
            mass_correction_term = 1 + (2 * (me_MeV / m_MeV) / gamma) + (me_MeV / m_MeV)^2;
            Wm_MeV = 2 * me_MeV * beta^2 / (1 - beta^2) / mass_correction_term;
            beta2 = beta^2;
            I_MeV = I_eV * 1e-6;
            m_to_cm = 100;
            MeV_to_J = 1.60217646e-13;
            I_J = I_MeV * MeV_to_J;

            % Stopping power
            restricted = false;
            if E_restricted_keV > 0.0 && (E_restricted_keV / 1000.0) < Wm_MeV
                restricted = true;
            end

            SN11 = (2.0 * me_MeV * beta2) / (1.0 - beta2);
            assert(I_MeV > 0);
            SN11 = SN11 / I_MeV;

            if restricted
                Wm_MeV = E_restricted_keV * 1e-3;
            end
            SN12 = Wm_MeV / I_MeV;

            SN2 = beta2;
            if restricted
                SN2 = SN2 / 2;
                SN2 = (SN2 + (1.0 - beta2) * Wm_MeV) / (4.0 * me_MeV);
            end

            delta = 0.0;
            if phase ~= phase_undefined
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

                plasma_energy_J = sqrt((4 * pi * electron_density_m3 * re) * Dirac_constant_J_s * c);

                C = 1.0 + 2.0 * log(I_J / plasma_energy_J);

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

            stopping_number = 0.5 * log(SN11 * SN12) - SN2 - SN3;

            if Z ~= 1
                effective_charge = Z * (1.0 - exp(-125.0 * beta / (Z^(2.0/3.0))));
            else
                effective_charge = 1.0 - exp(-125.0 * beta);
            end

            z = Z;

            energy_loss_leading_term_MeV_cm2_g = 0.307075 * (ZT / AT) * (z^2) / beta2;

            if E_MeV >= BETHE_LOWER_LIMIT_E_MEV_U
                energy_loss_MeV_cm = energy_loss_leading_term_MeV_cm2_g * stopping_number;
            else
               energy_loss_MeV_cm = 0;
            end

            energy_loss_keV_um = energy_loss_MeV_cm / 10; % [keV/μm] AND [MeV/mm]
            energy_loss_MeV_m = energy_loss_MeV_cm / 100; % [MeV/m]

            disp(['dE/dx is ', num2str(energy_loss_keV_um), ' [MeV/mm]']);
            disp(['current Energy is ', num2str(E_MeV), ' [MeV]']);
            disp(['---------------------------------------------------']);
        end
        
        function saveResults(obj)
            non_zero_mask = any(obj.positions ~= 0, 1);
            filtered_positions = obj.positions(:, non_zero_mask);
            csvwrite('trajectory.csv', filtered_positions');
            disp(['Non-zero positions have been saved!']);
            
            figure;
            plot(filtered_positions(1, :), filtered_positions(2, :), 'b.-', 'LineWidth', 1);
            title('Trajectory of Helium under Magnetic Field with Energy Loss');
            xlabel('X Position (cm)');
            ylabel('Y Position (cm)');
            grid on;
            axis equal;
        end
        
        function energy = getEnergyAtStep(obj, step)
            % Method to retrieve energy at a given step
            energy = obj.energies(step);
        end

        function displayStepRange(obj)
            % Method to display the number of steps with non-zero values
            non_zero_mask = any(obj.positions ~= 0, 1);
            num_non_zero_steps = sum(non_zero_mask);
            disp(['Number of steps with values: ', num2str(num_non_zero_steps)]);
        end
    end
end
