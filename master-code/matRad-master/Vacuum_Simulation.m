classdef Vacuum_Simulation
    properties
        % Constant variables
        mp_MeV = 938.272013;
        mp = 1.67262158e-27; % [kg]
        c = 299792458 * 100; % [mm/s]
        atomic_mass_unit_MeV_c2 = 931.494; % MeV/c^2
        q = 1.602e-19; % Elementary charge

        % Variables one might change
        E0_MeV; % Initial proton energy
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
        crossing_step; % Store the step where x crosses zero
    end
    
    methods
        function obj = Vacuum_Simulation(SPR_values, E0_MeV, B, grid_step, total_distance, initial_position, initial_velocity)
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
            obj.crossing_step = -1; % Initialize crossing step to -1 (not crossed yet)
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
            initial_x = obj.initial_position(1); % Store the initial x position
            
            for step = 2:obj.num_steps
                half_grid_step = obj.grid_step / 2;
                half_time_step = half_grid_step / norm(obj.velocities(1, step - 1));

                F_mid = - obj.q * cross(obj.velocities(:, step - 1), obj.B);
                a_mid = F_mid / obj.mp;

                mid_velocities = obj.velocities(:, step - 1) + a_mid * half_time_step;
                mid_positions = obj.positions(:, step - 1) + obj.velocities(:, step - 1) * half_time_step;

                energy_loss_keV_um = 1e-15;
                obj.energy_losses_adjusted(step) = energy_loss_keV_um;
                obj.E_MeV = obj.E_MeV - energy_loss_keV_um * obj.grid_step;
                obj.energies(step) = obj.E_MeV;

                gamma = 1 + obj.E_MeV / obj.atomic_mass_unit_MeV_c2;
                beta = sqrt(1 - 1/gamma^2);
                if gamma < 1
                    disp(['Gamma dropped below 1 at step ', num2str(step)]);
                    break;
                end

                F_end = - obj.q * cross(mid_velocities, obj.B);
                a_end = F_end / obj.mp;

                obj.velocities(:, step) = mid_velocities + a_end * half_time_step;
                obj.positions(:, step) = mid_positions + mid_velocities * half_time_step;
                obj.accelerations(:, step) = a_end;
                obj.times(step) = obj.times(step - 1) + 2 * half_time_step;

                if (obj.positions(1, step) - initial_x) * (obj.positions(1, step - 1) - initial_x) < 0
                    if abs(obj.positions(2, step) - obj.positions(2, 1)) > 1e-3
                        radius = abs(obj.positions(2, step)) / 2;
                        disp(['Proton crossed the initial x position at step ', num2str(step), ...
                              ' with a radius of ', num2str(radius), ' cm']);
                        obj.crossing_step = step; % Save the crossing step
                        break;
                    end
                end
            end
        end
        
        function saveResults(obj)
    non_zero_mask = any(obj.positions ~= 0, 1);
    filtered_positions = obj.positions(:, non_zero_mask);
    csvwrite('trajectory.csv', filtered_positions');
    disp(['Non-zero positions have been saved!']);
    
    figure;
    % Plot the proton trajectory
    plot(filtered_positions(1, :), filtered_positions(2, :), 'b.-', 'LineWidth', 1);
    hold on;
    
    % If the proton crossed the initial x position, plot a red dot
    if obj.crossing_step > 0
        plot(obj.positions(1, obj.crossing_step), obj.positions(2, obj.crossing_step), 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
        % Automatically place the legend in the best location
        legend('E = 50 MeV, B = 1.5 T', sprintf('Crossed x=0 at y=%.2f cm', obj.positions(2, obj.crossing_step)), 'Location', 'best', 'FontSize', 20);
    else
        % Automatically place the legend in the best location
        legend('E = 50 MeV, B = 1.5 T', 'Location', 'best');
    end
    
    title('Proton trajectory under magnetic field influence in a vacuum', 'FontSize', 20);
    xlabel('x [cm]', 'FontSize', 20);
    ylabel('y [cm]', 'FontSize', 20);
    grid on;
    axis equal;

    % Change font size of ticks
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 20; % Set font size for x-axis ticks
    ax.YAxis.FontSize = 20; % Set font size for y-axis ticks
end


        
        function energy = getEnergyAtStep(obj, step)
            energy = obj.energies(step);
        end

        function displayStepRange(obj)
            non_zero_mask = any(obj.positions ~= 0, 1);
            num_non_zero_steps = sum(non_zero_mask);
            disp(['Number of steps with values: ', num2str(num_non_zero_steps)]);
        end
    end
end