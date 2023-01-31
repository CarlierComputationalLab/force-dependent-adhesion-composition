% SensPipeline v3

%%
close all;
clear;

%% Specify directory to save data
directory = 'D:\Thesis project\Master Folder\Results\Data\';
date1 = char(datetime('today', 'format', 'ddMMyyyy'));
savedirectory = fullfile(directory,'Sensitivity analysis results',date1); 

%% Simulation duration 
t_start = 0; 
t_end = 600; %s
dt = 0.005; %s
tspan = [t_start:dt:t_end];
t_steps = numel(tspan);

%% Global constants

% Force exerted by myosin motors on actin filaments
conc_m = 4; % concentration of myosin motors in the cell (uM)
F_m = 2; %pN - Force exerted by a single myosin motor

% Force thresholds for slip bonds of the three complexes. 
F_th1 = 2; %pN
F_th2 = 2.5; %pN
F_th3 = 3; %pN

% Force threshold for vinculin binding
F_vb1 = 5; %pN
F_vb2 = 12; %pN

% Rate increase factor for Pcomp breakdown 
RIF_pcomp = 1; % 1 implies a 100 percent increase in the rate after signal crosses threshold

%% Spring constants

k_vinc = 0.25; %pN/nm

%% Signal threshold
sig_thresh = 0.1;

%% Functions to calculate force/time-dependent rate constants, fractional extension of clutch

k_unf_unloaded = 1.53690747659748;
k_unf_factor = 0.0501825624410996;
func_tal_unf = @(Force,threshold) ((Force<=threshold).*k_unf_unloaded.*exp(k_unf_factor.*(Force./threshold)) + (Force>threshold).*k_unf_unloaded.*exp(k_unf_factor));

% catch bond 
A = 2.51675514872035;
b = 0.106855913680857;
C = 0.000134394969149903;
d = 0.189677012719433;
func_catch = @(A,b,C,d,Force) A.*exp(-b.*Force) + C.*exp(d.*Force);

% Fractional extension of clutch when clutch+substrate is stretched 
func_extClutch = @(kclutch, ksub) (1-kclutch./(ksub+kclutch));

% SDRM factor
func_concRateConst = @(signalConcentration, signalThreshold, sigma) (signalConcentration>signalThreshold).*1 + (signalConcentration<=signalThreshold).*(sigma.*(signalConcentration));

%% Base rate consts
base_consts = zeros(36,1); % to store fixed rate constans 

% Pre-complex formation
base_consts(1) = 0.12; %k1f
base_consts(2) = 0.095; %k1r

% Seed formation 
base_consts(3) = 0.021; %k2f
base_consts(4) = 0.001*base_consts(3); %k2r

% Cluster formation 
base_consts(5) = 0.021; %+ 0.005*k_sub; %k3f
base_consts(6) = 0.001*base_consts(5); %k3r

% Seed reinforcement reverse rates
base_consts(10) = 0.001; %k7r
base_consts(11) = 0.001; %k8r

% Cluster reinforcement reverse rates
base_consts(15) = 0.001; %k12r
base_consts(16) = 0.001;%k13r

% Cluster formation from active seed
base_consts(18) = 0.001;%k14r
base_consts(20) = 0.001;%k15r
base_consts(22) = 0.001; %k16r

% Talin refolding - Inactive seeds
base_consts(24) = 0; %k17r 
base_consts(26) = 0; %k18r

% Talin refolding - Inactive clusters
base_consts(28) = 0; %k19r
base_consts(30) = 0; %k20r

% Cluster breakdown to seeds
base_consts(32) = 0; %k21r
base_consts(34) = 0; %k22r
base_consts(35) = 0.1; % k23_vmax 
base_consts(36) = 2.13; % k23_km

%% Initial concentrations

init_int = 1;
init_tal = 1;
init_vinc = 1;
init_sig = 1;

%% Substrate stiffness range
% k_sub_range = [0.01,0.1];
k_sub_range = [0.01, 0.02:0.02:10, 12.5, 15:5:30, 40, 50, 60, 80, 100, 150, 200, 250, 350, 500, 650, 800, 1000];

%% Sensitivity analysis
run('SensitivityAnalysisV3.m')

%% Saving results

mkdir(savedirectory);
save(strcat(savedirectory,'\SensResultsAndSensitivity.mat'),'sensitivity','sensData')

% mkdir('D:\Thesis project\Master Folder\Results\Data\Sensitivity analysis results\Sens-10012023')
% save('D:\Thesis project\Master Folder\Results\Data\Sensitivity analysis results\Sens-10012023\SensResultsAndSensitivity.mat','sensitivity','sensData')

% load('D:\Thesis project\Master Folder\Results\Data\Sensitivity analysis\OtherTrials\15052022-4\SensResults.mat')
% load('D:\Thesis project\Master Folder\Results\Data\Sensitivity analysis\20052022-includingSigThresh\SensResults.mat')

%% Unpacking struct
% k_sub_range = [0.01, 0.02:0.02:10, 12.5, 15:5:30, 40, 50, 60, 80, 100, 150, 200, 250, 350, 500, 650, 800, 1000];
% 
% vels_unpacked = NaN(numel(params), numel(ParamRange), numel(k_sub_range));
% 
% for i = 1:size(sensData,2)
%     vels_unpacked(i,:,:) = sensData(i).meanVels;
% end

