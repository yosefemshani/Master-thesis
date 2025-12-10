LETVector = full(dij.mLETDose{1,1}(:, 1));

% Dimensionen des Dosis-Gitters
x_dim = dij.doseGrid.dimensions(2); % Anzahl der x-Koordinaten
y_dim = dij.doseGrid.dimensions(1); % Anzahl der y-Koordinaten
z_dim = dij.doseGrid.dimensions(3); % Anzahl der z-Schichten

% Schritt 2: Wandle den 1D-Vektor in eine 3D-Matrix um
LETMatrix = reshape(LETVector, [y_dim, x_dim, z_dim]);

% Schritt 3: Beschränke auf die zweite z-Schicht
LETMatrix_z2 = LETMatrix(:,:,2); % Dosiswerte in der zweiten z-Schicht

% Schritt 4: Summiere die Dosiswerte entlang der y-Achse für jede x-Koordinate
depthLETCurve = sum(LETMatrix_z2, 1); % Summe über y für jede x-Koordinate

x = (1:length(depthLETCurve));

% yosef_x = load('yosef_x.mat');
% yosef_x = yosef_x.yosef_x;
% yosef_y = load('yosef_y.mat');
% yosef_y = yosef_y.yosef_y;

figure;
hold on;
plot(x*1.09375, depthLETCurve/max(depthLETCurve), 'b.-', 'LineWidth', 1);
plot(yosef_x, yosef_y/max(yosef_y), 'r.-', 'LineWidth', 1);
hold off;
title('Energy loss of proton beam');
xlabel('Range (mm)');
ylabel('Energy loss dE/dx (idk yet)');
grid on;
ylim([0,1]);