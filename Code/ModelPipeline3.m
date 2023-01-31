% Model 12 - for modifications

% Starts with clustering

% Precomplex formation
% 
% R1 =	(k1f)*int*tal*vinc - (k1r)*Pcomp
% 
% Seed formation
% R2 = 	(k2f)*Pcomp^25 - (k2r)*S1
% 
% Cluster formation
% R3 =	(k3f)*S1^2 - (k3r)*C1
% 
% Actin binding - seeds
% R4 = 	(k4f)*S1 - (k4r)*S1a 
% R5 =	(k5f)*S2 - (k5r)*S2a
% R6 =	(k6f)*S3 - (k6r)*S3a
% 
% Seed reinforcement
% R7 =	(k7f)*S1a*vinc^25 - (k7r)*S2a
% R8 =	(k8f)*S2a*vinc^25 - (k8r)*S3a
% 
% Actin binding - clusts
% R9 =	(k9f)*C1 - (k9r)*C1a
% R10 = (k10f)*C2 - (k10r)*C2a
% R11 =	(k11f)*C3 - (k11r)*C3a
% 
% Cluster reinforcement
% R12 = (k12f)*C1a*vinc^50 - (k12r)*C2a
% R13 =	(k13f)*C2a*vinc^50 - (k13r)*C3a
% 
% Cluster formation from active seed
% R14 = (k14f)*S1a^2 - (k14r)*C1a
% R15 =	(k15f)*S2a^2 - (k15r)*C2a
% R16 = (k16f)*S3a^2 - (k16r)*C3a
% 
% Talin refolding - Inactive Seed
% R17 =	(k17f)*S3 - (k17r)*S2*vinc^25
% R18 =	(k18f)*S2 - (k18r)*S1*vinc^25

% Talin refolding - Inactive cluster
% R19 =	(k19f)*C3 - (k19r)*C2*vinc^50
% R20 =	(k20f)*C2 - (k20r)*C1*vinc^50
% 
% Cluster breakdown to seed
% R21 =	(k21f)*C3 - (k21r)*S3^2
% R22 =	(k22f)*C2 - (k22r)*S2^2

% Signal molecule 
% 
% R23: Sig --> inSig, MM kinetics. d[inSig]/dt = V_max * [Sig]/(K_m + [Sig])
% d[Sig]/dt = V_max * [Sig]/(K_m + [Sig])

% Time/Force dependent rate constants
% k4r = var_rate(1);
% k5r = var_rate(2);
% k6r = var_rate(3);
% k7f = var_rate(4);
% k8f = var_rate(5);
% k9r = var_rate(6);
% k10r = var_rate(7);
% k11r = var_rate(8);
% k12f = var_rate(9);
% k13f = var_rate(10);

%  Fixed rate constants

% k1f = base_consts(1);
% k1r = base_consts(2);
% k2f = base_consts(3);
% k2r = base_consts(4);
% k3f = base_consts(5);
% k3r = base_consts(6);
% k4f = base_consts(7);
% k5f = base_consts(8);
% k6f = base_consts(9);
% k7r = base_consts(10);
% k8r = base_consts(11);
% k9f = base_consts(12);
% k10f = base_consts(13);
% k11f = base_consts(14);
% k12r = base_consts(15);
% k13r = base_consts(16);
% k14f = base_consts(17);
% k14r = base_consts(18);
% k15f = base_consts(19);
% k15r = base_consts(20);
% k16f = base_consts(21);
% k16r = base_consts(22);
% k17f = base_consts(23); var (1/f)?
% k17r = base_consts(24); 0
% k18f = base_consts(25); var (1/f)?
% k18r = base_consts(26); 0
% k19f = base_consts(27); var (1/f)?
% k19r = base_consts(28); 0
% k20f = base_consts(29); var (1/f)?
% k20r = base_consts(30); 0
% k21f = base_consts(31); var (1/f)? 
% k21r = base_consts(32); 0
% k22f = base_consts(33); var (1/f)? 
% k22r = base_consts(34); 0
% k23_vmax = base_consts(35); 
% k23_km = base_consts(36);

% % Concentrations
% int = concs(1);
% tal = concs(2);
% vinc = concs(3);
% Pcomp = concs(4);
% S1 = concs(5);
% S2 = concs(6);
% S3 = concs(7);
% S1a = concs(8);
% S2a = concs(9);
% S3a = concs(10);
% C1 = concs(11);
% C2 = concs(12);
% C3 = concs(13);
% C1a = concs(14);
% C2a = concs(15);
% C3a = concs(16);
% Sig = concs(17);
% inSig = concs(18);

%%
% close all;
clear;

condition = 'baseline'; % for no maturation set condition as 'NoMaturation'
date1 = char(datetime('today', 'format', 'ddMMyyyy'));
% date1 = '01122022';
savedirectory = fullfile('D:\Thesis project\Master Folder\Results\Figures',date1,'\',condition);
% mkdir(savedirectory);

%% Show/save graphs?

saveGraphs = 0; % should graphs be saved? 
showGraphs = 'on'; % should the figures be visible? on/off

%% NA assembly disassembly data
NA_assemblyData = readmatrix('D:/Thesis project/Master Folder/Calibration data/NA_assembly_disassembly_choi et al.csv');
NA_assemblyData_norm = NA_assemblyData;
NA_assemblyData_norm(:,2) = normalize(NA_assemblyData_norm(:,2),'range');

%% Simulation duration 
t_start = 0; 
t_end = 600; %s
dt = 0.005; %s
tspan = [t_start:dt:t_end];
t_steps = numel(tspan);

%% Global constants

% Concentration discretization 
cell_vol = 1; % cell vol in um^3 
factor_conc_to_molecules = cell_vol*602; % concentration in uM, 
    
% unloaded actin velocity 
v_u = 110; %nm/s

% Force exerted by myosin motors on actin filaments
conc_m = 4; % concentration of myosin motors in the cell (uM)
F_m = 2; %pN - Force exerted by a single myosin motor
F_myo = conc_m*factor_conc_to_molecules*F_m; %pN

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
k_tal = 0.1; %pN/nm
k_vinc = 0.25; %pN/nm
k_4tal = 4*k_tal; % stiffness of 1/4th talin molecule

% Effective spring constants of the complexes
kc1 = ((((k_4tal^-1 + k_4tal^-1 + k_4tal^-1)^-1) + k_vinc)^-1 + k_4tal^-1)^-1;
kc2 = (((((((k_4tal^-1 + k_4tal^-1)^-1) + k_vinc)^-1 + k_4tal^-1)^-1) + k_vinc)^-1 + k_4tal^-1)^-1;
kc3 = ((((((((k_4tal + k_vinc)^-1 + k_4tal^-1)^-1) + k_vinc)^-1 + k_4tal^-1)^-1) + k_vinc)^-1 + k_4tal^-1)^-1;

% Effective spring constants of the clusters
kseed1 = 25*kc1;
kseed2 = 25*kc2;
kseed3 = 25*kc3;
kclust1 = 2*kseed1;
kclust2 = 2*kseed2;
kclust3 = 2*kseed3;

k_clutches = [kseed1;kseed2;kseed3;kclust1;kclust2;kclust3];

% Stiffnesses of different sub-components of the clutch 
% m represents springs CDG, n represents springs BCDFG (refer supplementary
% material). 1,2, and 3 represent the components m and n for complex 1, 
% complex 2 and complex 3 respectively. 

k_m1 = k_4tal/2;
k_n1 = k_4tal/3;
k_m2 = k_4tal/2; 
k_n2 = ((1/k_4tal)+(1/(k_vinc + k_m2)))^-1;
k_m3 = ((1/k_4tal)+(1/(k_4tal + k_vinc)))^-1;
k_n3 = ((1/k_4tal)+(1/(k_vinc + k_m3)))^-1;

%% Max force each complex can experience

% Component of extension experienced by 4th talin subspring in each complex
% when an individual complex experiences an extension of X 

comp_c1 = (1-kc1/k_4tal)*(1-k_n1/k_4tal)*(1-k_m1/k_4tal); % *X
comp_c2 = (1-kc2/k_4tal)*(1-k_n2/k_4tal)*(1-k_m2/k_4tal);
comp_c3 = (1-kc3/k_4tal)*(1-k_n3/k_4tal)*(1-k_m3/k_4tal);

max_x4tal1 = F_th1/k_4tal; %nm  - Maximum stretch of 4th talin subspring. It reaches the slip bond force of 2pN at a stretch of 5nm if FL-talin stiffness is 0.1pN/nm. 
max_x4tal2 = F_th2/k_4tal;
max_x4tal3 = F_th3/k_4tal;

max_stretch_complex1 = max_x4tal1/comp_c1; %nm - Maximum stretch of the entire complex when the 4th talin subspring is at its maximum stretch
max_stretch_complex2 = max_x4tal2/comp_c2;
max_stretch_complex3 = max_x4tal3/comp_c3;

max_fc1 = kc1*max_stretch_complex1; % pN - Maximum force carrying capacity of the complex
max_fc2 = kc2*max_stretch_complex2;
max_fc3 = kc3*max_stretch_complex3;
max_fc = [max_fc1; max_fc2; max_fc3; max_fc1; max_fc2; max_fc3]; % storing the max force values

%% Signal threshold
sig_thresh = 0.1;

%% Functions to calculate force/time-dependent rate constants, fractional extension of clutch

% Talin unfolding
k_unf_unloaded = 1.53690747659748;
k_unf_factor = 0.0501825624410996;
func_tal_unf = @(Force,threshold) ((Force<=threshold).*k_unf_unloaded.*exp(k_unf_factor.*(Force./threshold)) + (Force>threshold).*k_unf_unloaded.*exp(k_unf_factor));

% Catch bond 
A = 2.51675514872035;
b = 0.106855913680857;
C = 0.000134394969149903;
d = 0.189677012719433;
func_catch = @(A,b,C,d,Force) A.*exp(-b.*Force) + C.*exp(d.*Force);% + 150.*exp(-Force./0.2);

% Slip bond
kslip_unloaded = 0.35; % unloaded slip bond rupture rate
func_kslip = @(Force,threshold) kslip_unloaded.*exp(Force./threshold);

% Fractional extension of clutch when clutch+substrate is stretched 
func_extClutch = @(kclutch, ksub) (1-kclutch./(ksub+kclutch));

% SDRM factor
func_concRateConst = @(signalConcentration, signalThreshold, sigma) (signalConcentration>signalThreshold).*1 + (signalConcentration<=signalThreshold).*(sigma.*(signalConcentration));

% SDRM factor if signal increase is considered instead of decay
func_concRateConst2 = @(signalConcentration, signalThreshold, sigma) (signalConcentration<signalThreshold).*1 + (signalConcentration>=signalThreshold).*(sigma.*(1-signalConcentration));

% TDRM parameters and function
k_sens = 0.05;
func_timeRateConst = @(time_bound) k_sens.*(time_bound)*dt;

% TDRM
TDRM_bool = 1; % if 1, TDRM is active. Set to 0 to make TDRM inactive

% Catch bond TDRM - Should it also be modified? 1 = yes; 0 = no 
cat = 0; 

%% Base rate consts
base_consts = zeros(36,1); % to store fixed rate constans 

% Assigning values to named variables (to be able to change them through a for loop)
talin_refold = 1;
talin_refold_factor = 0.5;
k_act = 1.5;

% Pre-complex formation
base_consts(1) = 0.12; %k1f
base_consts(2) = 0.095; %k1r

% Seed formation 
base_consts(3) = 0.021; %k2f
base_consts(4) = 0.001*base_consts(3); %k2r

% Cluster formation 
base_consts(5) = 0.021; %k3f
base_consts(6) = 0.001*base_consts(5); %k3r

% Actin binding rate - seeds
base_consts(7) = k_act; %k4f
base_consts(8) = k_act; %k5f
base_consts(9) = k_act; %k6f

% Seed reinforcement reverse rates
base_consts(10) = 0.001; %k7r
base_consts(11) = 0.001; %k8r

% Actin binding rate - clusters
base_consts(12) = k_act; %k9f
base_consts(13) = k_act; %k10f
base_consts(14) = k_act; %k11f

% Cluster reinforcement reverse rates
base_consts(15) = 0.001; %k12r
base_consts(16) = 0.001;%k13r

% Cluster formation from active seed
base_consts(17) = 1; %k14f
base_consts(18) = 0.001;%k14r
base_consts(19) = 1; %k15f
base_consts(20) = 0.001;%k15r
base_consts(21) = 1; %k16f
base_consts(22) = 0.001; %k16r

% Talin refolding - Inactive seeds
base_consts(23) = talin_refold_factor*talin_refold; %k17f
base_consts(24) = 0; %k17r 
base_consts(25) = talin_refold; %k18f
base_consts(26) = 0; %k18r

% Talin refolding - Inactive clusters
base_consts(27) = talin_refold_factor*talin_refold; %k19f
base_consts(28) = 0; %k19r
base_consts(29) = talin_refold; %k20f
base_consts(30) = 0; %k20r

% Cluster breakdown to seeds
base_consts(31) = 0.005; %k21f
base_consts(32) = 0; %k21r
base_consts(33) = 0.008; %k22f
base_consts(34) = 0; %k22r
base_consts(35) = 0.1; % k23_vmax 
base_consts(36) = 2.13; % k23_km

%% Initial concentrations

init_int = 1;
init_tal = 1;
init_vinc = 1;
init_sig = 1;

%% Substrate stiffness range
k_sub_range = [0.1,1,10,100];

%% Master storage variables

results_concs = NaN(length(k_sub_range),18, t_steps); 
results_vels = NaN(length(k_sub_range), t_steps);
results_forces = NaN(length(k_sub_range),7, t_steps);
results_varRates = NaN(length(k_sub_range), 10, t_steps);
results_FperSpecies = NaN(length(k_sub_range), 6, t_steps);
results_sig_dep_rates = NaN(length(k_sub_range), 12, t_steps);

%% for loop - Substrate stiffness

for k = 1:length(k_sub_range)
    
    disp([num2str(k), 'of', num2str(length(k_sub_range))])
    storeIter = 1;
    
    k_sub = k_sub_range(k);

    fractionalExt = func_extClutch(k_clutches,k_sub);
    
    %% Temporary storage variables
    
    temp_concs = zeros(size(results_concs,2),1);
%     temp_varRates = NaN(size(results_varRates,2),1);
    temp_forces = zeros(size(results_forces, 2),1);
%     temp_vels  = NaN(1,1);
%     temp_exts = NaN(6,1);
    
    % initializing variables to store values
    concs = NaN(size(results_concs,2),t_steps);% to store concentrations
    varRates = NaN(size(results_varRates,2),t_steps); % to store variable rate constants
    forces = NaN(size(results_forces,2),t_steps); % to store forces
    vels = NaN(1,t_steps); % to store retro velocities
    FperSpecies = NaN(6,t_steps);

    results_sig_dep_rates(k,:,1) = base_consts([1,2,3,4,5,6,7,12,25,29,31,33]);
    
    %% Variable to store clutch actin bound time
    t_bound = zeros(6,1);
    
    %% Set initial concentrations  
 
    concs(1,1) = init_int; % integrin
    concs(2,1) = init_tal; % talin
    concs(3,1) = init_vinc; % vinculin
    concs(17,1) = init_sig; % signal
%     concs(19,1) = init_sig;
    
    temp_concs(1) = concs(1,1);
    temp_concs(2) = concs(2,1);
    temp_concs(3) = concs(3,1);
    temp_concs(17) = concs(17,1);
%     temp_concs(19) = concs(19,1);
    
%% for loop - Integration

    consts = base_consts;

    for t = 2:length(tspan)
        
        % % SDRM
        factor = func_concRateConst(temp_concs(17),sig_thresh,1/sig_thresh);

        % Actin-binding of S1 and C1 reduces to 0 
        consts([7,12]) = base_consts(7)*factor;
        
        % Breakdown of S2 and C2 to S1 and C1 reduces to 0 
        consts([25,29]) = base_consts(25)*(factor);

        % Breakdown of AUB clusts to AUB seeds reduces to 0 
        consts([31,33]) =  base_consts([31, 33])*factor;

        % Breakdown of Pcomp increases by a factor of RIF_pcomp 
        consts([2]) = base_consts([2]) + (temp_concs(17)<=sig_thresh)*RIF_pcomp*base_consts([2])*(1-factor);
        
        % S1 and C1 disassembly rate is 0.0001*assembly rate (assembly rate = 0.021) in the beginning, increases to 2*0.01217 (choi et al.) 
        consts([4,6]) = (temp_concs(17)>sig_thresh)*base_consts([4,6]) + (temp_concs(17)<=sig_thresh)*2*0.01217*(1-factor);

        % Seed and cluster formation also reduce
        consts([1,3,5]) = base_consts([1,3,5])*factor;

        % Storing modified signal dependent rates
        results_sig_dep_rates(k,1:12,t) = consts([1,2,3,4,5,6,7,12,25,29,31,33]);
        
        slipUL = (TDRM_bool==0) + (TDRM_bool==1)*(1+func_timeRateConst(t_bound));
        
        % Force on individual complexes in the seeds and clusters

        F_ind = temp_forces(1:6)./[25; 25; 25; 50; 50; 50];
        
        % Calculating rate constants
        catchSlip_rates = func_kslip(F_ind,max_fc).*slipUL + func_catch(A.*((cat==0)+(cat==1).*slipUL),b,C.*((cat==0)+(cat==1).*slipUL),d,F_ind);
        talunf_rates = func_tal_unf(F_ind([1,2,4,5]),[F_vb1; F_vb2; F_vb1; F_vb2]); 
       
        % Storing force-dependent rates
        temp_varRates([1:3,6:8]) = catchSlip_rates;
        temp_varRates([4,5,9,10]) = talunf_rates;
     
        % Euler integration

        concs_new = temp_concs + DE_definitions11(temp_concs, temp_varRates, consts)*dt;
        concs_new(3) = concs(3,1); % Comment out if vinculin concentration should not be constant. 
        
        % Checking catch-slip bonds
        broken = find(F_ind>=max_fc);
        brokenSeeds = broken(broken<=3);
        brokenClusts = broken(broken>3);

        % Move AB clutch concs to AUB clutch concs
        concs_new(brokenSeeds+4) = concs_new(brokenSeeds+4)+concs_new(brokenSeeds+7);
        concs_new(brokenClusts+7) = concs_new(brokenClusts+7)+concs_new(brokenClusts+10);
        concs_new(brokenSeeds+7) = 0;
        concs_new(brokenClusts+10) = 0;
 
        % Reset force on recently broken AB clutches to 0
        temp_forces(broken) = 0;
        
        % Increase t_bound by 1
        t_bound = t_bound + 1;
        t_bound(broken) = 0;

        % Updating concentrations
        temp_concs = concs_new;
              
        % Calculate total force : number of molecules*force per molecule
        active_species_count = factor_conc_to_molecules*temp_concs([8,9,10,14,15,16]);
        force_per_species = active_species_count.*temp_forces(1:6);
        total_force = sum(force_per_species);

        % Store total force
        temp_forces(7) = total_force;
        
        % Calculate new retrograde velocity and displacement 
        v_retro = v_u*(1-(total_force/F_myo)); 
        x_disp4 = v_retro*dt;

        % Update force on each complex - reset to 0 when threshold reached, increase otherwise
        temp_forces(1:6) = (temp_concs([8:10, 14:16])>0).*(temp_forces(1:6)<=([25;25;25;50;50;50].*max_fc)).*(temp_forces(1:6) + fractionalExt.*(k_clutches.*x_disp4));

        forces(:,storeIter) = temp_forces;
        concs(:,storeIter) = concs_new;
        vels(:,storeIter) = v_retro;
        varRates(:,storeIter) = temp_varRates;
        FperSpecies(:,storeIter) = force_per_species;
        storeIter = storeIter + 1;
        
    end

    results_concs(k,:,:) = concs(:,:);
    results_vels(k,:) = vels(:);
    results_forces(k,:,:) = forces(:,:);
    results_varRates(k,:,:) = varRates(:,:);
    results_FperSpecies(k,:,:) = FperSpecies(:,:);
    
end

%% Res species
res_species = zeros(length(k_sub_range),16,t_steps);
for i = 1:length(k_sub_range)
    temp_concs_res = squeeze(results_concs(i,:,:));
    res_species(i,1,:) = 25*sum(temp_concs_res([5,6,7],:)); % AUB seed
    res_species(i,2,:) = 25*sum(temp_concs_res([8,9,10],:)); % AB seed
    res_species(i,3,:) = 50*sum(temp_concs_res([11,12,13],:)); % AUB clusts
    res_species(i,4,:) = 50*sum(temp_concs_res([14,15,16],:)); % AB clusts
    res_species(i,5,:) = 25*sum(temp_concs_res([6,7],:)); % AUB mid and higher order seeds
    res_species(i,6,:) = 50*sum(temp_concs_res([12,13],:)); % AUB mid and higher order clusters
    res_species(i,7,:) = 25*sum(temp_concs_res([9,10],:)); % AB mid and higher order seeds
    res_species(i,8,:) = 50*sum(temp_concs_res([15,16],:)); % AB mid and higher order clusters
    res_species(i,9,:) = sum(res_species(i,[5:8],:)); % all mid and high order species
    res_species(i,10,:) = 50*sum(temp_concs_res([13,16],:)); % AB and AUB C3
    res_species(i,11,:) = 50*sum(temp_concs_res([12,15],:)); % AB and AUB C2
    res_species(i,12,:) = 25*sum(temp_concs_res([7,10],:)); % AB and AUB S3
    res_species(i,13,:) = 25*sum(temp_concs_res([6,9],:)); % AB and AUB S2
    res_species(i,14,:) = sum(res_species(i,[10,12],:)); % AB and AUB high order species
    res_species(i,15,:) = 25*temp_concs_res(10, :) + 50*temp_concs_res(16,:); % AB S3 and C3
    res_species(i,16,:) = 25*temp_concs_res(9, :) + 50*temp_concs_res(15,:); % AB S2 and C2
end

%% Frequency of oscillations

% % Finding the time points where force hits 0
period = NaN(length(k_sub_range), 6, 85000);
zero_mat = NaN(size(period,2),size(period,3));

for k = 1:length(k_sub_range)
    forces = squeeze(results_forces(k,[1:6],:)); 
    for i = 1:size(forces,1)
        zero_idx = find(forces(i,1:end-1)<=0);
        zero_mat(i,1:length(zero_idx)) = zero_idx;
    end
    zero_sub = [zeros(6,1),zero_mat(:,1:end-1)];
    period(k,:,:) = (zero_mat-zero_sub)*dt;
end

period(period==dt) = NaN; % removing initial entries that are consecutively 0 for some species - these occur when the species are yet to be formed, so for the first few time steps. 
frequency = 1./(period); 

meanFreq = NaN(length(k_sub_range),6);
meanPeriod = NaN(length(k_sub_range),6);

for k = 1:length(k_sub_range)
    meanFreq(k,:) = mean(squeeze(frequency(k,:,:)),2,'omitnan')';
    meanPeriod(k,:) = mean(squeeze(period(k,:,:)),2,'omitnan')';
end


S1a = NaN(length(k_sub_range),1);
S2a = NaN(length(k_sub_range),1);
S3a = NaN(length(k_sub_range),1);
C1a = NaN(length(k_sub_range),1);
C2a = NaN(length(k_sub_range),1);
C3a = NaN(length(k_sub_range),1);

for k = 1:length(k_sub_range) 
    SubstrateStiffness{k,1} = num2str(k_sub_range(k));
    S1a(k,1) = meanPeriod(k,1);
    S2a(k,1) = meanPeriod(k,2);
    S3a(k,1) = meanPeriod(k,3);
    C1a(k,1) = meanPeriod(k,4);
    C2a(k,1) = meanPeriod(k,5);
    C3a(k,1) = meanPeriod(k,6);
end

T = table(SubstrateStiffness,S1a,S2a,S3a,C1a,C2a,C3a);

% % Time of concentration peak
% peakConcTime = NaN(length(k_sub_range),6);
% peakConc = NaN(length(k_sub_range),6);
% 
% rows = [8,9,10,14,15,16];
% for k = 1:length(k_sub_range)
%     for i = 1:length(rows)
%         p = rows(i);
%         [peakConc(k,i), peakConcTime(k,i)] = max(squeeze(results_concs(k,p,1:210/dt)));
%     end
% end
% peakConcTime = peakConcTime*dt;

%% Mass conservation check

nIDs = 8;
alphabet = ('a':'z').';
% alphabet = ["i", 'ii', 'iii', 'iv', 'v', 'vi', 'vii', 'viii'].';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('\bf(',chars,')'); 

% Mass cons integrin

f1h = figure;
sfh1 = subplot(1,3,[1]);
hold on

for k = 1:length(k_sub_range)
    concs = squeeze(results_concs(k,:,:));
    concs([5,6,7,8,9,10],:) = 25*concs([5,6,7,8,9,10],:);
    concs([11,12,13,14,15,16],:) = 50*concs([11,12,13,14,15,16],:);
    concs(:,end) = concs(:,end-1); % Without this the last value remains nan and leaves space in the plot
    plot(tspan, sum(concs([1,4,5:16],:)),'LineWidth',1.2,'DisplayName',strcat('k_s_u_b = ',num2str(k_sub_range(k))));
    data.IntMassCons(:,k) = sum(concs([1,4,5:16],:))';
end
subplotLabel(1) = text(0.025,0.95,charlbl{1},'Units','normalized','FontSize',14);
ylim([1-5e-14, 1+5e-14])
set(gca,'XGrid','on','YGrid','on','YTickLabel',...
    {'1 - 1e-13','1 - 5e-14','1','1 + 5e-14','1 + 1e-13'},...
    'XTick',[0,300,600])
title('Integrins', 'fontsize', 14)
xlabel('Time (s)', 'fontsize', 14)
ylabel('Concentration ($\mu$M)', 'fontsize', 14)

% % % % Mass cons Talin

sfh2 = subplot(1,3,[2]);
hold on
for k = 1:length(k_sub_range)
    concs = squeeze(results_concs(k,:,:));
    concs([5,6,7,8,9,10],:) = 25*concs([5,6,7,8,9,10],:);
    concs([11,12,13,14,15,16],:) = 50*concs([11,12,13,14,15,16],:);
    concs(:,end) = concs(:,end-1); % Without this the last value remains nan and leaves space in the plot
    plot(tspan, sum(concs([2,4,5:16],:)),'LineWidth',1.2,'DisplayName',strcat('k_s_u_b = ',num2str(k_sub_range(k))));
    data.TalMassCons(:,k) = sum(concs([2,4,5:16],:))';
end
subplotLabel(2) = text(0.025,0.95,charlbl{2},'Units','normalized','FontSize',14);
ylim([1-5e-14, 1+5e-14])
set(gca,'XGrid','on','YGrid','on','YTickLabel',[],...
    'XTick',[0,300,600])
title('Talin', 'fontsize', 14)
xlabel('Time (s)', 'fontsize', 14)

sfh4 = subplot(1,3,[3]);
plot(tspan, ones(4,length(tspan))*nan, 'LineWidth', 1.2)
set(gca, 'xtick', [],'visible','off')
leg = legend('$k{_s}{_u}{_b}$ = 0.1','$k{_s}{_u}{_b}$ = 1','$k{_s}{_u}{_b}$ = 10','$k{_s}{_u}{_b}$ = 100', 'Location', 'west');
title(leg,{"Substrate" + newline + "Stiffness" + newline + "(pN/nm)"})
% set(findall(subplotLabel, '-property', 'FontWeight'),'FontWeight','bold')
nicePlot(f1h, 14, 20,0.65, savedirectory, 'IntTalMassCons',saveGraphs, showGraphs)

% % % % % % Mass cons Vinculin
% % % % figure
% % % % sfh3 = subplot(1,3,3);
% % % % hold on
% % % % for k = 1:length(k_sub_range)
% % % %     concs = squeeze(results_concs(k,:,:));
% % % %     concs([5,8],:) = 1*25*concs([5,8],:);
% % % %     concs([6,9],:) = 2*25*concs([6,9],:);
% % % %     concs([7,10],:) = 3*25*concs([7,10],:);
% % % %     concs([11,14],:) = 1*50*concs([11,14],:);
% % % %     concs([12,15],:) = 2*50*concs([12,15],:);
% % % %     concs([13,16],:) = 3*50*concs([13,16],:);
% % % %     concs(:,end) = concs(:,end-1);
% % % %     plot(tspan, sum(concs([3,4,5:16],:)),'LineWidth',1.2,'DisplayName',strcat('k_s_u_b = ',num2str(k_sub_range(k))));
% % % % end
% % % % ylim([1-5e-14, 1+5e-14])
% % % % set(gca,'XGrid','on','YGrid','on','YTickLabel',...
% % % %     {'1 - 1e-13','1 - 5e-14','1','1 + 5e-14','1 + 1e-13'},...
% % % %     'XTick',[0,300,600])
% % % % title('Vinculin')% - mass conservation')
% % % % xlabel('Time (s)')
% % % % ylabel('Concentration (uM)')

% % % Plot signal decay and change of rate constant based on signal concentration 
concs = squeeze(results_concs(1,:,:));
sig_thresh_time_idx = (find(concs(17,:)<=sig_thresh,1));
sig_thresh_time = tspan(sig_thresh_time_idx);
f3h = figure;
hold on
plot(tspan, concs(17,:),'DisplayName','Signal concentration','LineWidth',1.2) 
plot(tspan, func_concRateConst(concs(17,:),0.1,1/0.1),'DisplayName', 'SDRM factor','LineWidth',1.2, 'color', '#EDB120');
% plot(tspan(sig_thresh_time_idx:end), func_concRateConst(concs(17,sig_thresh_time_idx:end),0.1,1/0.1),'DisplayName', 'SDRM factor','LineWidth',1.2, 'color', '#EDB120');
plot([0 sig_thresh_time] ,[concs(17,sig_thresh_time_idx) concs(17,sig_thresh_time_idx)],'r--','LineWidth',1.2,'DisplayName','Signal${_t}{_h}{_r}{_e}{_s}{_h}$');
plot([sig_thresh_time sig_thresh_time] ,[concs(17,sig_thresh_time_idx) 0],'DisplayName','$t{_s}{_i}{_g}$','LineWidth',1.2,'color','r');
legend();
% title('Signal concentration vs time')
xlabel('Time (s)', 'fontsize', 14);
ylabel('Concentration ($\mu$M)', 'fontsize', 14, 'interpreter', 'latex');
nicePlot(f3h, 14, 20, 0.65, savedirectory, 'SignalDecay',saveGraphs, showGraphs)


% % % Plot signal growth and change of rate constant based on signal
concs = squeeze(results_concs(1,:,:));
sig_thresh_time_idx = (find(concs(18,:)>=(1-sig_thresh),1));
sig_thresh_time = tspan(sig_thresh_time_idx);
f3h = figure;
hold on
plot(tspan, concs(18,:),'DisplayName','Signal concentration','LineWidth',1.2) 
plot(tspan, func_concRateConst2(concs(18,:),1-0.1,1/0.1),'DisplayName', 'SDRM factor','LineWidth',1.2, 'color', '#EDB120');
% plot(tspan(sig_thresh_time_idx:end), func_concRateConst(concs(17,sig_thresh_time_idx:end),0.1,1/0.1),'DisplayName', 'SDRM factor','LineWidth',1.2, 'color', '#EDB120');
plot([0 sig_thresh_time] ,[concs(18,sig_thresh_time_idx) concs(18,sig_thresh_time_idx)],'r--','LineWidth',1.2,'DisplayName','Signal${_t}{_h}{_r}{_e}{_s}{_h}$');
plot([sig_thresh_time sig_thresh_time] ,[concs(18,sig_thresh_time_idx) 0],'DisplayName','$t{_s}{_i}{_g}$','LineWidth',1.2,'color','r');
legend();
% title('Signal concentration vs time')
xlabel('Time (s)');
ylabel('Concentration ($\mu$M)');

nicePlot(f3h, 12, 19, 0.75, savedirectory, 'SignalGrowth',saveGraphs, showGraphs)

%% Concentration vs time - S1a, S2a, S3a, C1a, C2a, C3a for all stiffnesses

LineColor = {'0.1','1','10','100';'#0072BD', '#D95319', '#EDB120', '#7E2F8E'};
speciesFigsNames = ["S1", "S2", "S3", "S1a","S2a","S3a","C1", "C2", "C3", "C1a","C2a","C3a"];
Labels = ["(i)","(ii)","(iii)","(iv)","(v)","(vi)"];

tstart = 1;
tend = 600.005/dt;

% Arranging substrate stiffness from highest to lowest concentration of each species
clear idxMat
for i = 1:12
    concs = squeeze(results_concs(:,i+4,tstart:tend));
    [~,idxMat(i,:)] = sort(max(concs,[],2),'descend');
end

% save('D:/Thesis project/Master Folder/Code/idxMatrix.mat','idxMat');
% load('D:/Thesis project/Master Folder/Code/idxMatrix.mat','idxMat');

% % % %  Tiled plots for entire duration for all stiffnesses

plots = [5,6,7,11,12,13]+3;

f4h = figure;
t1 = tiledlayout(2,3);%tiledlayout(2,2);
for i = 1:length(plots)
    p = plots(i);
    stiffnessOrder = idxMat(p-4,:);
    speciesFigs{i} = nexttile;
    hold on
    for j = 1:length(stiffnessOrder)
        k = stiffnessOrder(j);
        concs = squeeze(results_concs(k,p,tstart:tend))';
        if p>4 && p<11
            concs = 25*concs;
        else 
            concs = 50*concs;
        end
        linecolor = LineColor{2,find(strcmp(LineColor(1,:),num2str(k_sub_range(k))))};
        plot(tspan(tstart:tend), concs,'LineWidth',1.2,'color',linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(k))," pN/nm"));
        xlim([tspan(tstart) tspan(tend)]);
    end
    title(strcat(speciesFigsNames(p-4),{' '}))
    subplotLabel(i) = text(0.4,0.95,Labels{i},'Units','normalized','FontSize',12);
end
leg = legend(speciesFigs{2}, 'Location', 'northoutside','orientation','horizontal','NumColumns',4);
% title(leg,"Substrate Stiffness (pN/nm)")
% title(leg,{"Substrate" + newline + "Stiffness" + newline + "(pN/nm)"})
% title(t1,'Concentration of integrins in AB species');%mid and high order species')
xlabel(t1,'Time (s)','Interpreter', 'latex', 'fontsize', 14)
ylabel(t1,'Concentration ($\mu$M)','Interpreter', 'latex', 'fontsize', 14)

figname = ['IntegrinsABspecies_allTime_',condition];
nicePlot(f4h, 12,  19, 0.75, savedirectory, figname,saveGraphs, showGraphs)

% % % %  Tiled plots for first 70 seconds for all stiffnesses
tend2 = 70/dt;
f4h = figure;
t1 = tiledlayout(2,3);
for i = 1:length(plots)
    p = plots(i);
    stiffnessOrder = idxMat(p-4,:);
    speciesFigs{i} = nexttile;
    hold on
    for j = 1:length(stiffnessOrder)
        k = stiffnessOrder(j);
        concs = squeeze(results_concs(k,p,tstart:tend2))';
        if p>4 && p<11
            concs = 25*concs;
        else 
            concs = 50*concs;
        end
        linecolor = LineColor{2,find(strcmp(LineColor(1,:),num2str(k_sub_range(k))))};
        plot(tspan(tstart:tend2), concs,'LineWidth',1.2,'color',linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(k))," pN/nm"));
        xlim([tspan(tstart) tspan(tend2)]);
    end
    title(strcat(speciesFigsNames(p-4),{' '}))
    subplotLabel(i) = text(0.4,0.95,Labels{i},'Units','normalized','FontSize',12);
end
leg = legend(speciesFigs{2}, 'Location', 'northoutside','orientation','horizontal','NumColumns',4);
% title(leg,"Substrate Stiffness (pN/nm)")
xlabel(t1,'Time (s)','Interpreter', 'latex', 'fontsize', 12)
ylabel(t1,'Concentration ($\mu$M)','Interpreter', 'latex', 'fontsize', 12)

figname = ['IntegrinsABspecies_first70',condition];
nicePlot(f4h, 12,  19, 0.75, savedirectory, figname, saveGraphs, showGraphs)

% % Maturation fraction (concentration of mid and high order seeds+clusts) over time
f1h = figure;
for i = 1:length(k_sub_range)
    k = k_sub_range(i);
    hold on
    concs = squeeze(results_concs(i,:,:));
    concs(5:10,tstart:tend)=concs(5:10,tstart:tend)*25;
    concs(11:16,tstart:tend)=concs(11:16,tstart:tend)*50;
    concs = squeeze(sum(concs([6,7,9,10,12,13,15,16],tstart:tend)));
    linecolor = LineColor{2,find(strcmp(LineColor(1,:),num2str(k_sub_range(i))))};
    plot(tspan(tstart:tend), concs,'LineWidth',1.2,'color', linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(i))))
end
leg = legend('Location', 'eastoutside');
title(leg,{"Substrate" + newline + "Stiffness" + newline + "(pN/nm)"})
xlabel('Time (s)','Interpreter', 'latex', 'fontsize', 14)
ylabel('Concentration ($\mu$M)','Interpreter', 'latex', 'fontsize', 14)

figname = ['MF_',condition];
nicePlot(f1h, 14, 19, 0.75, savedirectory, figname, saveGraphs, showGraphs)

% % Concentration of low, mid and high order seeds+clusts over time
f2h = figure;
t1 = tiledlayout(1,3);
fig_low = nexttile(t1,1);
ylim([0, 0.5]);
title('Low order')
fig_mid = nexttile(t1,2);
ylim([0, 0.5]);
title('Mid order')
fig_high = nexttile(t1,3);
ylim([0, 0.5]);
title('High order')

for i = 1:length(k_sub_range)
    k = k_sub_range(i);
    concs = squeeze(results_concs(i,:,:));
    concs(5:10,tstart:tend)=concs(5:10,tstart:tend)*25;
    concs(11:16,tstart:tend)=concs(11:16,tstart:tend)*50;
    concs_low = squeeze(sum(concs([5,8,11,14],tstart:tend)));
    concs_mid = squeeze(sum(concs([6,9,12,15],tstart:tend)));
    concs_high = squeeze(sum(concs([7,10,13,16],tstart:tend)));
    linecolor = LineColor{2,find(strcmp(LineColor(1,:),num2str(k_sub_range(i))))};
    hold(fig_low,'on')
    plot(fig_low, tspan(tstart:tend), concs_low,'LineWidth',1.2,'color', linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(i))))

    hold(fig_mid,'on')
    plot(fig_mid, tspan(tstart:tend), concs_mid,'LineWidth',1.2,'color', linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(i))))

    hold(fig_high,'on')
    plot(fig_high, tspan(tstart:tend), concs_high,'LineWidth',1.2,'color', linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(i))))
end
leg = legend('Location', 'eastoutside');
title(leg,{"Substrate" + newline + "Stiffness" + newline + "(pN/nm)"})
xlabel(t1, 'Time (s)','Interpreter', 'latex', 'fontsize', 14)
ylabel(t1,'Concentration ($\mu$M)','Interpreter', 'latex', 'fontsize', 14)

figname = ['ClutchConcByOrder_',condition];
nicePlot(f2h, 14, 22, 0.45, savedirectory, figname, saveGraphs, showGraphs)

% % Combined plot of the above two figures
f2h = figure;
t1 = tiledlayout(2,3);
fig_low = nexttile(t1,1);
ylim([0, 0.5]);
title('Low order')
fig_mid = nexttile(t1,2);
ylim([0, 0.5]);
title('Mid order')
fig_high = nexttile(t1,3);
ylim([0, 0.5]);
title('High order')
fig_MF = nexttile(t1,[1,3]);
title('Mature adhesions')
% fig_2_1 = nexttile(t1,4);
% fig_2_1.Visible = 'off';
% 
% fig_2_3 = nexttile(t1,6);
% fig_2_3.Visible = 'off';


for i = 1:length(k_sub_range)
    k = k_sub_range(i);
    concs = squeeze(results_concs(i,:,:));
    concs(5:10,tstart:tend)=concs(5:10,tstart:tend)*25;
    concs(11:16,tstart:tend)=concs(11:16,tstart:tend)*50;
    concs_low = squeeze(sum(concs([5,8,11,14],tstart:tend)));
    concs_mid = squeeze(sum(concs([6,9,12,15],tstart:tend)));
    concs_high = squeeze(sum(concs([7,10,13,16],tstart:tend)));
    concs_MF = squeeze(sum(concs([6,7,9,10,12,13,15,16],tstart:tend)));
    linecolor = LineColor{2,find(strcmp(LineColor(1,:),num2str(k_sub_range(i))))};
    hold(fig_low,'on')
    plot(fig_low, tspan(tstart:tend), concs_low,'LineWidth',1.2,'color', linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(i))," pN/nm"))
    
    hold(fig_mid,'on')
    plot(fig_mid, tspan(tstart:tend), concs_mid,'LineWidth',1.2,'color', linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(i))," pN/nm"))
        
    hold(fig_high,'on')
    plot(fig_high, tspan(tstart:tend), concs_high,'LineWidth',1.2,'color', linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(i))," pN/nm"))
        
    hold(fig_MF,'on')
    plot(fig_MF,tspan(tstart:tend), concs_MF,'LineWidth',1.2,'color', linecolor,'DisplayName',strcat('$k{_s}{_u}{_b}$ = ',num2str(k_sub_range(i))," pN/nm"))
end
leg = legend(fig_mid, 'Location', 'northoutside','orientation','horizontal','NumColumns',4);
% title(leg,"Substrate Stiffness (pN/nm)")
xlabel(t1, 'Time (s)','Interpreter', 'latex', 'fontsize', 12)
ylabel(t1,'Concentration ($\mu$M)','Interpreter', 'latex', 'fontsize', 12)

figname = ['MF+ClutchConcByOrder_',condition];
nicePlot(f2h, 12, 19, 0.75, savedirectory, figname, saveGraphs, showGraphs)

% %

% % Concentration of S1 and C1 when no maturation takes place
if strcmp(condition,'NoMaturation')
    concs = squeeze(results_concs(1,:,:));
    sig_thresh_time_idx = (find(concs(17,:)<=sig_thresh,1));
    sig_thresh_time = tspan(sig_thresh_time_idx);

    % % % Concentration of seeds and clusts over time with no maturation
    f2h = figure;
    hold on
    for k = 1%:length(k_sub_range)
        concs = squeeze(results_concs(k,:,:));
        concs(5,:) = 25*concs(5,:);
        concs(11,:) = 50*concs(11,:); 
        plot(tspan, sum(concs([5,11],:)), 'DisplayName', '[S1]+[C1]' , 'LineWidth',1.2);%strcat('k_s_u_b = ',num2str(k_sub_range(k)))
    end
    concBelow0point1_time_idx = (find(sum(concs([5,11],:))>=0.1,1,'last'));
    concAbove0point1_time_idx = (find(sum(concs([5,11],:))>=0.1,1));
    concBelow0point1_time = tspan(concBelow0point1_time_idx);
    concAbove0point1_time = tspan(concAbove0point1_time_idx);
    abovetime = concBelow0point1_time - concAbove0point1_time;
    % title('Concentration of integrins in seeds and clusts over time with no maturation')
    xlabel('Time (s)', 'fontsize', 14);
    ylabel('Concentration ($\mu$M)', 'fontsize', 14);
    plot(NA_assemblyData_norm(:,1),max(sum(concs([5,11],:)))*NA_assemblyData_norm(:,2),'DisplayName', 'Experimental data (Choi et al. (2008))')
    plot([sig_thresh_time, sig_thresh_time], [0,sum(concs([5,11],sig_thresh_time_idx))],'--','color','r','DisplayName','$t{_s}{_i}{_g}$'); 
    plot([concAbove0point1_time concAbove0point1_time], [0, sum(concs([5,11],concAbove0point1_time_idx))], '--', 'color', 'k', 'HandleVisibility', 'off');
    plot([concBelow0point1_time concBelow0point1_time], [0, sum(concs([5,11],concBelow0point1_time_idx))], '--', 'color', 'k', 'HandleVisibility', 'off');
    an1 = annotation(f2h,'doublearrow',[0.1375 0.341785714285714],[0.16852380952381 0.16852380952381]);
    an2 = annotation('textbox',[0.2 0.168523809523813 0.204285714285714 0],'EdgeColor','none','String',sprintf('%0.1f s',abovetime),'FontSize',13);
    legend
    
    figname = ['SeedsClusts',condition];
    nicePlot(f2h,14,20,0.65, savedirectory, figname,saveGraphs, showGraphs);
end

%% Force exerted by each species
speciesNames = ["S1a","S2a","S3a","C1a","C2a","C3a"];

% Fitting a curve over the peaks of the force-time curves for each species
FperSp_smooth = NaN(length(k_sub_range), 6, length(tspan));
total_force_smooth = NaN(length(k_sub_range), length(tspan));

for k = 1:length(k_sub_range)
    FperSpecies = squeeze(results_FperSpecies(k, :, :));
    for i = 1:size(FperSpecies,1)
        np = ceil(meanPeriod(k,i)/dt);
        [a,b] = findpeaks(FperSpecies(i,:),'MinPeakDistance',0.8*np); % a is the y value, b is the time point
        p = fit([0,tspan(b)]',[0,a]','smoothingspline', 'SmoothingParam',0.002);
        FperSp_smooth(k,i,:) = feval(p,tspan);
        FperSp_smooth(FperSp_smooth <0) = 0;
    end
    [c,d] = findpeaks(sum(FperSpecies), 'MinPeakDistance',5000);
    q = fit([0, tspan(d)]',[0,c]','smoothingspline','SmoothingParam',0.0002);
    total_force_smooth(k,:) = feval(q,tspan);
end

% % % Tiled version: 
f3h = figure;
t3 = tiledlayout(1,4);
f3h.Position = f3h.Position.*[0.1,1,2,1];
% t3.Position = t3.Position.*[1,1,1.1,1];
for k = 1:length(k_sub_range)
    FperSpecies = squeeze(results_FperSpecies(k, :, :));
    FperSpec_plot{k} = nexttile; 
    hold on
    ylim([0 (ceil(max(total_force_smooth(:))/100)*100)+200])
    for i = 1:size(FperSpecies,1)
        plot(tspan, squeeze(FperSp_smooth(k,i,:))', 'DisplayName', speciesNames(i),'LineWidth',1.2) 
        title(strcat('$k{_s}{_u}{_b}$ = ',{' '}, num2str(k_sub_range(k)),{' '}, 'pN/nm'));
    end
    plot(tspan, total_force_smooth(k,:) ,'--', 'DisplayName','Total force','LineWidth',1.3)    
    subplotLabel(k) = text(0.4,0.95,Labels{k},'Units','normalized','FontSize',12);
end
legend('Location','eastoutside')
t3.Padding = 'compact';
% title(t3, strcat('Force exerted by each species on different stiffnesses'));
xlabel(t3, 'Time (s)','Interpreter', 'latex', 'fontsize', 12);
ylabel(t3, 'Force (pN)','Interpreter', 'latex', 'fontsize', 12);

figname = ['ForceBySpecies_',condition];
nicePlot(f3h, 12, 19, 0.40, savedirectory, figname, saveGraphs, showGraphs)

% % % Plot the same separately for each stiffness. 
% for k = 1:length(k_sub_range)
%     FperSpecies = squeeze(results_FperSpecies(k, :, :));
%     figure 
%     hold on
%     for i = 1:size(FperSpecies,1)
%         plot(tspan, squeeze(FperSp_smooth(k,i,:))', 'DisplayName', speciesNames(i),'LineWidth',1.2) 
%         title(strcat('ksub = ', num2str(k_sub_range(k))));
%         legend
%     end
%     hold on
%     plot(tspan, total_force_smooth(k,:) ,'--', 'DisplayName','Total force','LineWidth',1.3)    
%     legend
%     title(strcat('Force exerted by each species for k_s_u_b = ',num2str(k_sub_range(k))));
%     xlabel('Time (s)');
%     ylabel('Force (pN)');
%     ylim([0 ceil(max(total_force_smooth(:))/1000)*1000])
% end

%% Plot rate-constants over time

% % Time dependent rate constants
rateNames = ["$k{_4}{_r}$", "$k{_5}{_r}$", "$k{_6}{_r}$", "$k{_7}{_f}$", "$k{_8}{_f}$", "$k{_9}{_r}$", "$k{_1}{_0}{_r}$", "$k{_1}{_1}{_r}$", "$k{_1}{_2}{_f}$", "$k{_1}{_3}{_f}$"];

% % K4r and k11r were chosen as they correspond to the catch-slip bond 
% % rupture rate of the sfotest (S1a) and stiffest species (C3a) respectively
ratePlots = [ "$k{_4}{_r}$","$k{_1}{_1}{_r}$" ];

stiffness = [0.1,100]; % list of stiffnesses for which plots are required
plot_time_lim = [200, 5]; % time scale of plot: 0 to plot_time_lim, specified in the same order as stiffnesses

f6h = figure;
t = tiledlayout(2,2);
c = 1;
for i = 1:length(ratePlots)
    p = find(rateNames == ratePlots(i));
    for j = 1:length(stiffness)
        k = find(k_sub_range == stiffness(j));
        plotData = squeeze(results_varRates(k, p, :));
        sfh{c} = nexttile;
        plot(tspan(1:plot_time_lim(j)/dt), plotData(1:plot_time_lim(j)/dt), 'DisplayName', rateNames(p))
        hold on
        plot([tspan(1) tspan(plot_time_lim(j)/dt)], [mean(plotData(1:plot_time_lim(j)/dt)) mean(plotData(1:plot_time_lim(j)/dt))], 'DisplayName', 'Mean', 'LineWidth', 1);
        xlim([0 plot_time_lim(j)])
        if c <=2
            title(strcat('$k{_s}{_u}{_b}$ = ', {' '}, num2str(stiffness(j)),{' '}, 'pN/nm'));
        end
        subplotLabel(c) = text(0.025,0.95,charlbl{c},'Units','normalized','FontSize',12);
        c=c+1;
    end
end
ylabel(t,'Rate (1/s)', 'Interpreter', 'latex')
xlabel(t, 'Time (s)', 'Interpreter', 'latex')
% title(t, 'With TDRM', 'Interpreter', 'latex')

figname = ['ReactionRates_',condition];
nicePlot(f6h, 14, 20, 0.65, savedirectory, figname, saveGraphs, showGraphs) 

% % Signal dependent rate constants 
% colors = linspecer(9);
% rateNames = ["k1f", "k1r", "k2f", "k2r","k3f","k3r","k4f", "k9f", "k18f", "k20f","k21f","k22f"];
% ratePlots = [1,2,3,4,5,6,7,12,25,29,31,33];
%  
% plots = [1,3,5,7,12,25,29,31,33];[7];
% 
% figure
% hold on
% for i = 1:length(plots)
%     p = plots(i);
%     q = find(ratePlots == p);
%     plot(tspan, squeeze(results_sig_dep_rates(k,q,:)),'DisplayName',rateNames(q),'LineWidth',1.2);%, 'color', colors(i,:));
% end
% title('Signal-dependent rate decrease');
% xlabel('Time (s)')
% ylabel('Rate (1/s)')
% l=legend('Location', 'east')
%% Plot mean velocity for each stiffness
meanV = mean(results_vels(1:4,:),2,'omitnan');
figure
plot(k_sub_range, meanV)
set(gca,'xscale','log')
