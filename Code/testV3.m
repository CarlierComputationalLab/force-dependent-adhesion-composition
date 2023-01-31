% 'test' version 3

% file that runs the simulation for a given set of conditions
% Uses for

% parameters that change :
% 'cell_vol', 'v_u', 'k_tal', 'kslip_unloaded', 'k_sens', ...
% 'init_int', 'init_tal', 'init_vinc', 'init_sig', ...
% 'F_th1', 'F_th2', 'F_th3', ...
% 'RIF_pcomp', ...
% 'k14f', 'k15f', 'k16f', 'k21f', 'k22f', 'talin_refold', 'talin_refold_factor', 'k_act', 'sig_thresh'


%%

factor_conc_to_molecules = cell_vol*602; % concentration in uM, cell vol in um^-3
F_myo = conc_m*factor_conc_to_molecules*F_m; %pN
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

% stiffnesses of different sub-components of the clutch
k_m1 = k_4tal/2;
k_n1 = k_4tal/3;
k_m2 = k_4tal/2; 
k_n2 = ((1/k_4tal)+(1/(k_vinc + k_m2)))^-1;
k_m3 = ((1/k_4tal)+(1/(k_4tal + k_vinc)))^-1;
k_n3 = ((1/k_4tal)+(1/(k_vinc + k_m3)))^-1;

%% Max force each complex can experience

comp_c1 = (1-kc1/k_4tal)*(1-k_n1/k_4tal)*(1-k_m1/k_4tal); % *X
comp_c2 = (1-kc2/k_4tal)*(1-k_n2/k_4tal)*(1-k_m2/k_4tal);
comp_c3 = (1-kc3/k_4tal)*(1-k_n3/k_4tal)*(1-k_m3/k_4tal);

max_x4tal1 = F_th1/k_4tal; %nm  - Maximum stretch of talin subspring that binds actin. It reaches the slip bond force of 2pN at a stretch of 5nm. 
max_x4tal2 = F_th2/k_4tal;
max_x4tal3 = F_th3/k_4tal;

max_stretch_complex1 = max_x4tal1/comp_c1;
max_stretch_complex2 = max_x4tal2/comp_c2;
max_stretch_complex3 = max_x4tal3/comp_c3;

max_fc1 = kc1*max_stretch_complex1;
max_fc2 = kc2*max_stretch_complex2;
max_fc3 = kc3*max_stretch_complex3;
max_fc = [max_fc1; max_fc2; max_fc3; max_fc1; max_fc2; max_fc3]; % storing the max force values

% maximum stretch of seed and cluster is the same as the corresponding
% complexes (because springs connected in parallel experience same
% extension)

%% Functions to calculate force/time-dependent rate constants
func_kslip = @(Force,threshold) kslip_unloaded.*exp(Force./threshold);
func_timeRateConst = @(time_bound) k_sens.*(time_bound)*dt;

%% Base constants which are checked in sensitivity analysis

% Actin binding rate - seeds
base_consts(7) = k_act; %k4f
base_consts(8) = k_act; %k5f
base_consts(9) = k_act; %k6f

% Actin binding rate - clusters
base_consts(12) = k_act; %k9f
base_consts(13) = k_act; %k10f
base_consts(14) = k_act; %k11f

% Cluster formation from active seed
base_consts(17) = k14f; %k14f
base_consts(19) = k15f; %k15f
base_consts(21) = k16f; %k16f

% Talin refolding - Inactive seeds
base_consts(23) = talin_refold_factor*talin_refold; %k17f
base_consts(25) = talin_refold; %k18f

% Talin refolding - Inactive Clusts
base_consts(27) = talin_refold_factor*talin_refold; %k19f
base_consts(29) = talin_refold; %k20f

% Cluster breakdown to seeds
base_consts(31) = k21f; %k21f
base_consts(33) = k22f; %k22f

%% Defining master storage vars

results_concs = NaN(length(k_sub_range),18,t_steps);
results_vels = NaN(length(k_sub_range),t_steps);
results_mForces = NaN(length(k_sub_range),t_steps);

%% for loop - Substrate stiffness

parfor k = 1:length(k_sub_range)
    
    disp([num2str(k), 'of', num2str(length(k_sub_range))])
    storeIter = 1;
    
    k_sub = k_sub_range(k);

    fractionalExt = func_extClutch(k_clutches,k_sub);

    %% Initialising storage variables
    
%     temp_concs = zeros(size(results_concs,2),1);
    temp_concs = zeros(18,1);
    temp_forces = zeros(7,1);
    temp_vels  = NaN(1,1);
    
    % initializing variables to store values
%     concs = NaN(size(results_concs,2),t_steps);% to store concentrations
    concs = NaN(18,t_steps);% to store concentrations
    vels = NaN(1,t_steps); % to store retro velocities
    mForces = NaN(1,t_steps);
    
    %% Variable to store clutch actin bound time
    t_bound = zeros(6,1);
    
    %% Set initial concentrations  
 
    concs(1,1) = init_int; % integrin
    concs(2,1) = init_tal; % talin
    concs(3,1) = init_vinc; % vinculin
    concs(17,1) = init_sig; % signal
    
    temp_concs(1) = concs(1,1);
    temp_concs(2) = concs(2,1);
    temp_concs(3) = concs(3,1);
    temp_concs(17) = concs(17,1);
    
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
 
        % % Time-dependent rate constant modification
       
        % % TDRM
        TDRM_bool = 1; % if 1, TDRM is active. Set to 0 to make TDRM inactive
        t_thresh = 0;
        
        slipUL = (TDRM_bool==0) + (TDRM_bool==1)*(1+func_timeRateConst(t_bound));

        % Force on individual complexes in the seeds and clusters

        F_ind = temp_forces(1:6)./[25; 25; 25; 50; 50; 50];
        
        % Catch bond TDRM - Should it also be modified? 1 = yes; 0 = no 
        cat = 0; 
        
        % Calculating rate constants
        catchSlip_rates = func_kslip(F_ind,max_fc).*slipUL + func_catch(A.*((cat==0)+(cat==1).*slipUL),b,C.*((cat==0)+(cat==1).*slipUL),d,F_ind);
        talunf_rates = func_tal_unf(F_ind([1,2,4,5]),[F_vb1; F_vb2; F_vb1; F_vb2]); 
       
        % Storing force-dependent rates
%         temp_varRates = [catchSlip_rates(1:3)', talunf_rates(1:2)', catchSlip_rates(4:6)',talunf_rates(3:4)']
        temp_varRates = NaN(10,1);
        temp_varRates([1:3,6:8]) = catchSlip_rates;
        temp_varRates([4,5,9,10]) = talunf_rates;
       
        % Euler integration

%         concs_new = temp_concs +  DE_definitions11_parfor(temp_concs, consts, k4r, k5r, k6r, k7f, k8f, k9r, k10r, k11r, k12f, k13f)*dt; 
        concs_new = temp_concs + DE_definitions11(temp_concs, temp_varRates, consts)*dt;
        concs_new(3) = concs(3,1);
        
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

        concs(:,storeIter) = concs_new;
        vels(:,storeIter) = v_retro;
        mForces(:,storeIter) = temp_forces(7); % store the total force in each timestep
        storeIter = storeIter + 1;

    end

    results_concs(k,:,:) = concs(:,:);
    results_vels(k,:) = vels(:);
    results_mForces(k,:) = mForces(:); % store total force in each timestep for each stiffness

end

%% find optimum stiffness for given parameter values
mean_vels_j = mean(results_vels,2,'omitnan')'; %find the average over time of velocity for each stiffness  - outputs vector of dimensions (numel(k_sub_range),1)
mean_forces_j = mean(results_mForces, 2, 'omitnan')'; %find the average over time of total force for each stiffness - outputs vector of dimensions (numel(k_sub_range),1)
[min_vels_j, minVelIdx_j] = min(mean_vels_j);
% [minVelStiffness,~] = ind2sub(size(results_vels),minVelIdx);
OptStiff_temp_j = k_sub_range(minVelIdx_j);


%% Find integrins in different species for given parameter values

res_species = zeros(length(k_sub_range),14,t_steps);
for i = 1:length(k_sub_range)
    temp_concs_res = squeeze(results_concs(i,:,:));
    res_species(i,1,:) = 25*sum(temp_concs_res([5,6,7],:)); % inactive seed
    res_species(i,2,:) = 25*sum(temp_concs_res([8,9,10],:)); % active seed
    res_species(i,3,:) = 50*sum(temp_concs_res([11,12,13],:)); % inactive cluster
    res_species(i,4,:) = 50*sum(temp_concs_res([14,15,16],:)); % active cluster
    res_species(i,5,:) = 25*sum(temp_concs_res([6,7],:)); % inactive mid and higher order seeds
    res_species(i,6,:) = 50*sum(temp_concs_res([12,13],:)); % inactive mid and higher order clusters
    res_species(i,7,:) = 25*sum(temp_concs_res([9,10],:)); % active mid and higher order seeds
    res_species(i,8,:) = 50*sum(temp_concs_res([15,16],:)); % active mid and higher order clusters
    res_species(i,9,:) = sum(res_species(i,[5:8],:)); % all mid and high order species
    res_species(i,10,:) = 50*sum(temp_concs_res([13,16],:)); % AB and AUB C3
    res_species(i,11,:) = 50*sum(temp_concs_res([12,15],:)); % AB and AUB C2
    res_species(i,12,:) = 25*sum(temp_concs_res([7,10],:)); % AB and AUB S3
    res_species(i,13,:) = 25*sum(temp_concs_res([6,9],:)); % AB and AUB S2
    res_species(i,14,:) = sum(res_species(i,[10,12],:)); % AB and AUB high order species
    res_species(i,15,:) = 25*temp_concs_res(10, :) + 50*temp_concs_res(16,:); % AB S3 and C3
    res_species(i,16,:) = 25*temp_concs_res(9, :) + 50*temp_concs_res(15,:); % AB S2 and C2
end

% Integrins in mid, high species - outputs vector of dimensions (numel(k_sub_range),1)
IntMHS_temp_j = res_species(:,9, end-1); 
