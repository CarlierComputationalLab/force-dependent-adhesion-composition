% SensitivityAnalysis version 3

params= {'cell_vol', 'v_u', 'k_tal', 'kslip_unloaded', 'k_sens', ...
    'init_int', 'init_tal', 'init_vinc', 'init_sig', ...
    'F_th1', 'F_th2', 'F_th3', ...
    'RIF_pcomp', ...
    'k14f', 'k15f', 'k16f', 'k21f', 'k22f', 'talin_refold', 'talin_refold_factor', 'k_act', 'sig_thresh'};

base_values = [1, 110, 0.1, 0.35, 0.05,...
    1, 1, 1, 1, ...
    2, 2.5, 3, ...
    1, ...
    1, 1, 1, 0.005, 0.008, 1, 0.5, 1.5, 0.1];

ParamRange = [0.8,0.9,1.1,1.2];

%%
sensData = struct([]); % the data is filled in such that each column is one stiffness, each row is a ParamRange value. Size of the struct is 1xnumel(params). There are 4 fields.  

%%

for ii =1:numel(params)
%%
    disp(['Parameter number ', num2str(ii), 'of', num2str(length(params))])
    restoreDefault(params,base_values)
    param = params{ii};
    
    IntMHS_temp = NaN(numel(ParamRange),numel(k_sub_range));
    mean_forces = NaN(numel(ParamRange),numel(k_sub_range));
    mean_vels = NaN(numel(ParamRange),numel(k_sub_range));
    OptStiff_temp = NaN(numel(ParamRange),1);
%%
    for jj = 1:numel(ParamRange)
        value = base_values(ii)*ParamRange(jj);
        assignin('base',param,value);
        run('testV3.m')
        IntMHS_temp(jj,:) = IntMHS_temp_j;
        mean_vels(jj,:) = mean_vels_j; 
        mean_forces(jj,:) = mean_forces_j;
        OptStiff_temp(jj) = OptStiff_temp_j;
        
        clear results_concs results_vels results_mForces res_species minVelStiffness  
    end
    sensData(ii).IntMHS = IntMHS_temp;
    sensData(ii).meanVels = mean_vels;
    sensData(ii).meanForces = mean_forces;
    sensData(ii).OptStiff = OptStiff_temp;
end

%% Running baseline simulation

restoreDefault(params, base_values);
run('testV3.m')

%% storing basline data in the sensitivity results

for ii = 1:numel(params)
    sensData(ii).IntMHS(end+1,:) = IntMHS_temp_j;
    sensData(ii).meanVels(end+1,:) = mean_vels_j;
    sensData(ii).meanForces(end+1,:) = mean_forces_j;
    sensData(ii).OptStiff(end+1,:) = OptStiff_temp_j;
end


%% Analysis

sensitivity = struct([]);
for ii = 1:numel(params)
    IntMHSdata = sensData(ii).IntMHS;
    OptStiffdata = sensData(ii).OptStiff;
        
    sens_IntMHS = (abs(IntMHSdata(1:end-1,:)-IntMHSdata(end,:))./(IntMHSdata(end,:)))./(abs(1-ParamRange(1:end))');
    sens_OptStiff = (abs(OptStiffdata(1:end-1)-OptStiffdata(end))./(OptStiffdata(end)))./(abs(1-ParamRange(1:end))');
       
    sensitivity(ii).IntMHS = sens_IntMHS;
    sensitivity(ii).OptStiff = sens_OptStiff;
end


%%
function restoreDefault(params,base_values) 
    for i =1:numel(params)
        assignin('base',params{i},base_values(i))
    end
end

