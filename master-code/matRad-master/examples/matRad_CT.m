%% Irradiating given CT patient data set with raw source code instead of GUI

%% Patient Data Import

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

load('TEST123.mat'); % TEST123 is DoubleSlice but with Prostate as CTV with Objective 50 Gy Mean Dose [defined in cst]

%% Treatment Plan

pln.radiationMode = 'protons';        
pln.machine       = 'Generic';
pln.propOpt.bioOptimization = 'const_RBExD';
pln.propDoseCalc.calcLET = 1;

pln.numOfFractions        = 30;
pln.propStf.gantryAngles  = [90 270];
pln.propStf.couchAngles   = [0 0];
pln.propStf.bixelWidth    = 1;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% Lets generate dosimetric information by pre-computing dose influence 
% matrices for unit beamlet intensities. Having dose influences available 
% allows for subsequent inverse optimization. 
dij = matRad_calcParticleDose(ct,stf,pln,cst);

%% Inverse Optimization for IMPT
% The goal of the fluence optimization is to find a set of bixel/spot 
% weights which yield the best possible dose distribution according to the 
% clinical objectives and constraints underlying the radiation treatment
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Let's plot the transversal iso-center dose slice
slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
figure
imagesc(resultGUI.RBExDose(:,:,slice)),colorbar,colormap(jet)

plane = 3;
doseWindow = [0 max([resultGUI.RBExDose(:);])];

figure,title('original plan')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.RBExDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);

%% Plot the Resulting Beam Dose Slice
% Let's plot the transversal iso-center dose slice of beam 1 and beam 2
% separately 
%figure
%subplot(121),imagesc(resultGUI.RBExDose_beam1(:,:,slice)),colorbar,colormap(jet),title('dose of beam 1')
%subplot(122),imagesc(resultGUI.RBExDose_beam2(:,:,slice)),colorbar,colormap(jet),title('dose of beam 2')

%%
% Now let's simulate a patient shift in y direction for both beams
%stf(1).isoCenter(2) = stf(1).isoCenter(2) - 4;
%stf(2).isoCenter(2) = stf(2).isoCenter(2) - 4;
%pln.propStf.isoCenter       = reshape([stf.isoCenter],[3 pln.propStf.numOfBeams])';

%% Recalculate Plan
% Let's use the existing optimized pencil beam weights and recalculate the RBE weighted dose
%resultGUI_isoShift = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);