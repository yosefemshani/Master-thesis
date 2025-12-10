%% Bethe Bloch formula for the calculation of energy loss of per distance (here: for protons) with the help of libamtrack
%% Author: Yosef Emshani, 02.02.2024 18:26

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
Z=1; %for protons

% variables one might change!
E_MeV_u = 70;
E_restricted_keV = 1000;
I_eV = 75; % 75 eV for liquid water. for others, see https://github.com/libamtrack/library/blob/44dc48cfa977c05008ad12646d798d6c4b6ea504/include/AT_DataMaterial.h l. 197-203
AT = 18.0006; % mass number for H2O
ZT = 10; % ordinary number for target: water (H2O: 2 * 1 + 1 * 8 = 10)
phase = phase_condensed; %for liquid water! for others, see https://github.com/libamtrack/library/blob/master/include/AT_DataMaterial.h#L246
                                                            %lines 240-253 for different materials!

%simple calculations and conversion factors
gamma = 1 + E_MeV_u / atomic_mass_unit_MeV_c2; %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L45 lines 65-82)
assert(gamma >= 1.0);

beta = sqrt(1 - 1/gamma^2); %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L45 line 75)
mass_correction_term = 1 + (2 * (me_MeV / m_MeV) / gamma) + (me_MeV / m_MeV)^2; %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L278 l. 281)
Wm_MeV = 2 * me_MeV * beta^2 / (1 - beta^2) / mass_correction_term; %(see https://github.com/libamtrack/library/blob/master/src/AT_PhysicsRoutines.c#L285 l.288)
beta2 = beta^2;
I_MeV	= I_eV * 1e-6;
m_to_cm = 100;
MeV_to_J = 1.60217646e-13;
I_J = I_MeV * MeV_to_J;

%% Stopping power

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
        gamma_single =  1.0 + E_MeV_u / atomic_mass_unit_MeV_c2;
        assert(gamma_single >= 1.0);
        
        assert(E_MeV_u >= 0);
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

    assert(SN11 > 0);
    assert(SN12 > 0);

    stopping_number = 0.5 * log(SN11 * SN12) - SN2 - SN3;

%% Leading energy loss term
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

%% Energy loss
% Compute only above 1.0 MeV, otherwise theory is too wrong below return zero
% TODO: Find smarter criterion because this may cause problems in the code (as it did
% TODO: with the inappropriately set lower limit for CSDA range integration (was 0, now 1.0 MeV)
if E_MeV_u >= BETHE_LOWER_LIMIT_E_MEV_U 
    energy_loss_MeV_cm = energy_loss_leading_term_MeV_cm2_g * stopping_number;
else
    energy_loss_MeV_cm = 0;
end

% unit conversion
energy_loss_keV_um = energy_loss_MeV_cm / 10;


% TEST OUTPUT DISPLAY FOR CANDIDATE VARIABLES
disp(['dE/dx is ', num2str(energy_loss_MeV_cm), ' [MeV/cm]']);
disp(['dE/dx is ', num2str(energy_loss_keV_um), ' [keV/μm]']);
