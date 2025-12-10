%% Irradiating given CT patient data set with raw source code instead of GUI

%% Patient Data Import

matRad_rc; %If this throws an error, run it from the parent directory first to set the paths

%load('TEST123.mat');
load('2_04_P_CT_CST.mat'); % DAS SELBE WIE TEST123, aber nur das CT!

% Ich will hier zum Beispiel
% C:\Users\MP_Bernoulli\Downloads\master-code\2_04_P_real_cropped_single\DICOM_DoubleSlice
% einlesen
% und der Prostata als einziges Organ als Target auswählen  STOP STOP falsche Directory !!!mit einer
% Objective: Mean Dose 50 Gy zum beispiel! danke habibi okay mach richtige
% hier bro was hast du alles für objectives reingeamcht hahaha der beam
% sieht anders aus als ich den kenne DU HAST NUR EIN FELD, damit wird
% schwer, zu optiieren, man braucht zwei 90 und 270 Grad HAST DU GESEHEN!!!
% meinst du das dvh komplett arsch ist? Ja aber normal, was erwartets du
% BESSER ne hahahah ja deutlich aber hiermit kriegen wir dann 200000000
% bixel spots bro dann rechne ich 200000 jahre auch. Deswegen hatte ich nur
% ein beam und nur eine objective prostate 50 gy mean dose. Okay und wo
% rechnest du genau? was genau "wo" und was? "Bro dann rechne ich 20000
% Jahrre"  hWoi ? hier: GARNICHT MEHR STEIL. Ja checke, aber ich frage mich
% inwiefern das für meine Arbeit wichtig ist, weißt du? Mir geht es ja
% nicht darum einen klinischn Sinnvollen Plann zu erstellen. Ja du
% brauchstt das alles net. Keine Ahnung was wirhier mache
% Der algorithmus braucht 2-3h für 88 bixel spots nur. Ja 2 Beams sind
% doppelt so lange hahah ja deswegen momentan noch schwer für development
% phase. wenn der code 100% funktioniert, schalte ich zweite beam an mit
% deinen objectives. Aber woher hast du die genauen zahlen für die Voobn
% jecCtihrvieSTian Doktorarbeit achs okay verstehe.ABE R ich habe selber
% variert, wichtig ist nur, dass die dosis erschränkungen ineghaoten werdne
% wir redenncoh mal drüber WRATE
% lass mich mal was ausprobieren. Ich habe 'TEST123' geladen und alles
% außer das CT und die Konturen also CST entfernt und das dann als neue
% .mat gespeichert Ge?Nau Ok tmm weil ich wollte iwann z.B. neues CT
% datensatz einlesen und dafür dachte ich könnte man hier eif CT ordner
% einlesen und er macht zu .mat weißt duCheke. Habe das bisher immer über
% GUI gemacht mit dem umwandeln von DICOM in .mat aber das geht safe warte
% kp erstaml kann man aber rausfinden easy ja scheiss drauf das reicht
% erstmal danke für deie hilfe brokein Ding wie lang bleibst noch büro?6
% wie 66Uhr digga was??Ja hahah bro... Ich mache mir jetzt Reis in der
% Mikrowelle dikka... guten hunger ich geh schlafen habibiGute Nacht
% bruder. Diesen Text musst du archivieren hahahah mach ich buderherz. Gute
% Nach Gutte nacht.

%% Objectives for OAR and TARGET

% for i = 1:size(cst, 1)
%     if strcmp(cst{i, 2}, 'Rectum')
%         cst{i, 3}                        = 'OAR';
%         cst{i, 5}.Priority               = 1;
%         cst{i, 6}{1, 1}.className        = 'DoseObjectives.matRad_MaxDVH';
%         cst{i, 6}{1, 1}.parameters{1, 1} = 2; 
%         cst{i, 6}{1, 1}.parameters{1, 2} = 60;
%         cst{i, 6}{1, 1}.penalty          = 1000;
%     end
% end

% % Keine Blase in den Slices vorhanden, scheinbar!
% for i = 1:size(cst, 1)
%     if strcmp(cst{i, 2}, 'Urinary bladder')
%         cst{i, 3}                        = 'OAR';
%         cst{i, 5}.Priority               = 1;
%         cst{i, 6}{1, 1}.className        = 'DoseObjectives.matRad_MaxDVH';
%         cst{i, 6}{1, 1}.parameters{1, 1} = 1; 
%         cst{i, 6}{1, 1}.parameters{1, 2} = 50;
%         cst{i, 6}{1, 1}.penalty          = 1000;
%     end
% end

% for i = 1:size(cst, 1)
%     if strcmp(cst{i, 2}, 'Femoral head L')
%         cst{i, 3}                        = 'OAR';
%         cst{i, 5}.Priority               = 2;
%         cst{i, 6}{1, 1}.className        = 'DoseObjectives.matRad_MaxDVH';
%         cst{i, 6}{1, 1}.parameters{1, 1} = 35; 
%         cst{i, 6}{1, 1}.parameters{1, 2} = 5;
%         cst{i, 6}{1, 1}.penalty          = 100;
%     end
% end

% for i = 1:size(cst, 1)
%     if strcmp(cst{i, 2}, 'Femoral head R')
%         cst{i, 3}                        = 'OAR';
%         cst{i, 5}.Priority               = 2;
%         cst{i, 6}{1, 1}.className        = 'DoseObjectives.matRad_MaxDVH';
%         cst{i, 6}{1, 1}.parameters{1, 1} = 35; 
%         cst{i, 6}{1, 1}.parameters{1, 2} = 5;
%         cst{i, 6}{1, 1}.penalty          = 100;
%     end
% end

% for i = 1:size(cst, 1)
%     if strcmp(cst{i, 2}, 'External')
%         cst{i, 3}                        = 'OAR';
%         cst{i, 5}.Priority               = 2;
%         cst{i, 6}{1, 1}.className        = 'DoseObjectives.matRad_MaxDVH';
%         cst{i, 6}{1, 1}.parameters{1, 1} = 61; % verschriebene Dosis in Gy
%         cst{i, 6}{1, 1}.parameters{1, 2} = 0;
%         cst{i, 6}{1, 1}.penalty          = 10000;
%     end
% end

for i = 1:size(cst, 1)
    if strcmp(cst{i, 2}, 'Prostate')
        cst{i, 3}                        = 'TARGET';
        cst{i, 5}.Priority               = 1; % Lass mich nochmal kurz schauen.
        cst{i, 6}{1, 1}.className        = 'DoseObjectives.matRad_MeanDose';
        cst{i, 6}{1, 1}.parameters{1,1}  = 50; % verschriebene Dosis in Gy
        cst{i, 6}{1, 1}.penalty          = 1; 
    end
end

% Das ist jetzt nur Target und External Maximum. Wollte nur sichergehen,
% wegen WORKSPACE CST. warum eigentlich rbexdose undnicht physical dose.
% Spielt eig keien Rolle. Wir nehmen einen constRBE von 1.1 an. Daher juckt
% nicht, aber ist ja jetzt auch kein Porblem. Stimmt. und die werte für
% external und prostata haben wir aus christian arbeit nh? Also die
% verschriebene Dosis sind 60 Gy. Wir wollen mindestens 95% unsere CTVs mit
% 95% von 60 Gy vollballern. So und die Kurve im DVH wollen wir möglichst
% steil hinbekommen, yani einen sigmoidalen Verlauf. daher wählen wir das
% External Maximum so niedrig wie möglich. Ich zeige dir mal was sonst
% passiert. WARTE
% das ist ja so maximal arsch
% das it ok ich würd das so lassen erstmal weil man da auch gut das CTV rot
% ausgemalt hat! ich kann dir jetzt spontan net sagen, womit du optimale
% sachen kriegst du musst aufjedenfall mit "PENALATY" und Dosis spielen
% Was genau sind denn die werte von christian 1.
% und 2. wie kann ich meandose 50 Gy wieder machen z.b:? Christans werte
% sind überflüassig für dich. Da gilt für die Riskoorgane. Aber ich mach
% dir 50 warte
% Löwe wirklich
% Ok letzte frage woher hast du das nochmal gemacht

% TEST123 is DoubleSlice but with Prostate as CTV with Objective 50 Gy Mean Dose [defined in cst]

%% Treatment Plan

pln.radiationMode = 'protons';        
pln.machine       = 'Generic';
pln.propOpt.bioOptimization = 'const_RBExD';
pln.propDoseCalc.calcLET = 1;

pln.numOfFractions        = 1; 
pln.propStf.gantryAngles  = [270];
pln.propStf.couchAngles   = [0];
pln.propStf.bixelWidth    = 5;
pln.propStf.numOfBeams    = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter     = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);

csvwrite('isoCenter.csv', pln.propStf.isoCenter');
disp(['Iso Center positions have been saved as isoCenter.csv!']);

%pln.propStf.isoCenter     = [0,0,0];
pln.propOpt.runDAO        = 0;
pln.propOpt.runSequencing = 0;

% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Current CTV PBS grid (targetPoint_bev) with iso center at ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% [-40,10000,-2], [-40,10000, 0] [-40,10000, 2] [-40,10000,4]
% [-38,10000,-2], [-38,10000, 0] [-38,10000, 2] [-38,10000,4]
% ...
% [38,10000,4]

%% TEST adjusting stf parameters (rayPos e.g.) 
% Apparently, changing rayPos, rayPos_bev and targetPoint changes NOTHING

%stf(1).ray(1).rayPos = [0, 50, 50];
%stf(1).ray(1).rayPos_bev = [50, 0, 50];

%for i = 1:length(stf(1).ray)
%    stf(1).ray(i).rayPos = [0, 0, 0];
%end

%for i = 1:length(stf(1).ray)
%    stf(1).ray(i).rayPos_bev = [0, 0, 0];
%end

%for i = 1:length(stf(1).ray)
%   stf(1).ray(i).targetPoint_bev = [-40, 10000, -2]; % [y,z(?)/SAD(?),x]
%end

%for i = 1:length(stf(1).ray)
%    stf(1).ray(i).rayPos_bev = [-20,0,-1];
%end

%% Extract targetPoint_bev and save to CSV
% Initialize a matrix to store targetPoint_bev values
%allTargetPoints = [];

% Loop through each ray in stf(1)
%for i = 1:length(stf(1).ray)
    % Extract targetPoint_bev for each ray
%    targetPoint = stf(1).ray(i).targetPoint_bev;

    % Append to the matrix (concatenate vertically)
%    allTargetPoints = [allTargetPoints; targetPoint];
%end

% Save the matrix to a CSV file named 'target.csv'
%writematrix(allTargetPoints, 'target.csv');

% Or use csvwrite (deprecated in newer versions of MATLAB, use writematrix instead)
% csvwrite('target.csv', allTargetPoints);


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