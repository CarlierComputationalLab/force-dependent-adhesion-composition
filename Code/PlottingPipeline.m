% Plotting V2

%% Save directory
condition = 'baseline';
date1 = char(datetime('today', 'format', 'ddMMyyyy'));
% date1 = '10012023';
savedirectory = fullfile('D:\Thesis project\Master Folder\Results\Figures',date1,condition,'sensFigs');
mkdir(savedirectory);

%% saveGraphs?
saveGraphs = 0; % should figures be saved? 0/1
showGraphs = 'on'; % should figures be visible? on/off

%% load data

load('D:\Thesis project\Master Folder\Results\Data\Sensitivity analysis results\Sens-13012023\SensResultsAndSensitivity.mat')

%% Params and stiffness
params= {'$cell-vol$', '$v{_u}$', 'k${_t}{_a}{_l}$', 'k${_s}{_l}{_i}{_p}{_{{_U}{_L}}}$', 'k${_s}{_e}{_n}{_s}$', ...
    'Initial${_i}{_n}{_t}$', 'Initial${_t}{_a}{_l}$', 'Initial${_v}{_i}{_n}{_c}$', 'Initial${_s}{_i}{_g}$', ...
    'F${_t}{_h}{_{_1}}$', 'F${_t}{_h}{_{_2}}$', 'F${_t}{_h}{_{_3}}$', ...
    'RIF${_{pcomp}}$', ...
    'k${_1}{_4}{_f}$', 'k${_1}{_5}{_f}$', 'k${_1}{_6}{_f}$', 'k${_2}{_1}{_f}$', 'k${_2}{_2}{_f}$', 'tal${_r}{_f}$', 'tal${_r}{_f}{_{{_f}{_a}{_c}{_t}{_o}{_r}}}$', 'k${_a}{_c}{_t}$', 'sig${_t}{_h}{_r}{_e}{_s}{_h}$'};

% Plots

% plotParams = {'$v{_u}$', 'k${_t}{_a}{_l}$', 'k${_s}{_l}{_i}{_p}{_{{_U}{_L}}}$', 'k${_s}{_e}{_n}{_s}$', ...
%     'Initial${_i}{_n}{_t}$', 'Initial${_v}{_i}{_n}{_c}$', ...
%     'F${_t}{_h}{_{_1}}$',...
%     'tal${_r}{_f}$', 'tal${_r}{_f}{_{{_f}{_a}{_c}{_t}{_o}{_r}}}$', 'k${_a}{_c}{_t}$'};

plotParams = {'$v{_u}$', 'k${_t}{_a}{_l}$', ...
    'Initial${_v}{_i}{_n}{_c}$','Initial${_i}{_n}{_t}$', ...
    'k${_s}{_e}{_n}{_s}$', 'k${_a}{_c}{_t}$', ...
    'tal${_r}{_f}$', 'tal${_r}{_f}{_{{_f}{_a}{_c}{_t}{_o}{_r}}}$'};%

qlist = []; %stores the position indices of the selected params (plotParams) in the variable 'params'
for p = 1:length(plotParams)
    IndexC = strfind(params, plotParams(p));
    q = find(not(cellfun('isempty',IndexC))); 
    qlist = [qlist,q];
end

% k_sub_range = [0.01,0.1,1,10,100];
k_sub_range = [0.01, 0.02:0.02:10, 12.5, 15:5:30, 40, 50, 60, 80, 100, 150, 200, 250, 350, 500, 650, 800, 1000];
k_sub_sens = [0.1, 1, 10, 100]; %stiffnesses at which sensitivity for MF to different parameters should be plotted in a bar plot
[k_sub_diffs, k_sub_sens_idx] = min(abs(repmat(k_sub_range', [1, length(k_sub_sens)]) - k_sub_sens));

ParamRange = [0.8,0.9,1.1,1.2,1];

% Colors to be used in plots
colors = linspecer(length(params));

%% Data Prep
nIDs = 8;
alphabet = ('a':'z').';
chars = num2cell(alphabet(1:nIDs));
chars = chars.';
charlbl = strcat('\bf(',chars,')'); 

% I have 9 parameters to plot. 
% number of rows of yData is number of groups. So i need 4 rows, one for each stiffness - so this will be a 4x9 matrix
% And i have 4 such yDatas, one for each ParamRange
for r = 1:4 % paramRange
    for p = 1:length(params) %cycle through parameters
        IndexC = strfind(params,params(p));
        q = find(not(cellfun('isempty',IndexC)));
        yData_OptStiff(p) = sensitivity(q).OptStiff(r);
        for k = 1:length(k_sub_sens_idx) % stiffness 
            kk = k_sub_sens_idx(k);
            yData_IntMHS(k,p) = sensitivity(q).IntMHS(r,kk); %  for a given range value and stiffness, extract sensitivities of all parameters of interest.     
        end
    end
    yData_IntMHS_byRange{r} = yData_IntMHS;
    yData_OptStiff_byRange{r} = yData_OptStiff;
end

%% Param sensitivity plot - mean120_IntMHS - all params

f1h = figure;
t = tiledlayout(4,1);
% Top bar graph
for r = 1:4
    yData = yData_IntMHS_byRange{r};
    sens_IntMHS_byRange(r) = nexttile;
    plt = bar(sens_IntMHS_byRange(r),yData,'FaceColor', 'flat');
    for ii = 1:length(params)
        plt(ii).CData = ones(4,1)*colors(ii,:);
    end
    subplotLabel(r) = text(0.025,0.85,charlbl{r},'Units','normalized','FontSize',12);

    xl = get(gca, 'xlim');
    hold on
    plot(gca,xl, [1 1], 'r', 'lineWidth', 1)
    if ((1-ParamRange(r))*100) > 0
        ylabel(strcat({'$+'}, num2str((1-ParamRange(r))*100), {' '}, '\%$'),'Interpreter','latex');
    else
        ylabel(strcat({'$'}, num2str((1-ParamRange(r))*100), {' '}, '\%$'),'Interpreter','latex');
    end
    sens_IntMHS_byRange(r).YLim = [0 3.2];
    grid(sens_IntMHS_byRange(r), 'on')
    if r ~= 4
        sens_IntMHS_byRange(r).XAxis.Visible = 'off';
    end
    sens_IntMHS_byRange(r).XTickLabel = {'0.1', '1', '10', '100'};
end
f1h.Position = f1h.Position.*[1,1,1.5,1];
% leg = legend(sens_IntMHS_byRange(2),params, 'Location','northeastoutside');
leg = legend(sens_IntMHS_byRange(1),params, 'Location','northoutside','Orientation','horizontal', 'Numcolumns',6);
% leg = legend(sens_IntMHS_byRange(1),params, 'Location','northeastoutside','Orientation','vertical', 'Numcolumns',1);
xlabel(t, 'Substrate stiffness (pN/nm)','Interpreter','latex', 'fontsize', 14)
ylabel(t, 'Parameter sensitivity','Interpreter','latex', 'fontsize', 14)
% title(t, 'Percentage of integrins in mid and high order species at equlibrium', 'Interpreter', 'latex')
figname = ['Sens_IntMHS_allParams_',condition];
nicePlot(f1h, 14, 20, 0.75, savedirectory, figname, saveGraphs, showGraphs)

%% Param sensitivity plot - mean120_IntMHS - select params

f1h = figure;
t = tiledlayout(4,1);
% Top bar graph
for r = 1:4 % cycle through parameter range values (0.8, 0.9, 1.1, 1.2)
    yData = yData_IntMHS_byRange{r}; % returns a 4x22 matrix --> 4 stiffnesses (rows), 22 parameters (colums). 
    yData = yData(:,qlist); % subsetting data for the selected parameters
    sens_IntMHS_byRange(r) = nexttile;
    plt = bar(sens_IntMHS_byRange(r),yData,'FaceColor', 'flat');
    for ii = 1:length(plotParams)
        plt(ii).CData = ones(4,1)*colors(ii,:);
    end
    subplotLabel(r) = text(0.025,0.85,charlbl{r},'Units','normalized','FontSize',12);

%     plt.CData(
    xl = get(gca, 'xlim');
    hold on
    plot(gca,xl, [1 1], 'r', 'lineWidth', 1)
    if ((1-ParamRange(r))*100) > 0
        ylabel(strcat({'$+'}, num2str((1-ParamRange(r))*100), {' '}, '\%$'),'Interpreter','latex');
    else
        ylabel(strcat({'$'}, num2str((1-ParamRange(r))*100), {' '}, '\%$'),'Interpreter','latex');
    end
    sens_IntMHS_byRange(r).YLim = [0 3.2];
    grid(sens_IntMHS_byRange(r), 'on')
    if r ~= 4
        sens_IntMHS_byRange(r).XAxis.Visible = 'off';
    end
    sens_IntMHS_byRange(r).XTickLabel = {'0.1', '1', '10', '100'};
end
f1h.Position = f1h.Position.*[1,1,1.5,1];
% t.Position = t.Position.*[1,1,1.1,1];
% leg = legend(sens_IntMHS_byRange(2),plotParams, 'Location','northeastoutside');
leg = legend(sens_IntMHS_byRange(1),plotParams, 'Location','northoutside','Orientation','horizontal', 'Numcolumns',6);
% leg = legend(sens_IntMHS_byRange(1),plotParams, 'Location','northeastoutside','Orientation','vertical', 'Numcolumns',1);
xlabel(t, 'Substrate stiffness (pN/nm)','Interpreter','latex', 'fontsize', 14)
ylabel(t, 'Parameter sensitivity','Interpreter','latex', 'fontsize', 14)
% title(t, 'Percentage of integrins in mid and high order species at equlibrium', 'Interpreter', 'latex')

figname = ['Sens_IntMHS_selectParams_',condition];
nicePlot(f1h, 14, 20, 0.75,  savedirectory, figname ,saveGraphs, showGraphs)

%% Param sensitivity plot - OptStiff - all Params

f2h = figure;
t = tiledlayout(4,1);
for r = 1:4
    yData = yData_OptStiff_byRange{r};
    sens_OptStiff_byRange(r) = nexttile;
    hold on
    for n = 1:numel(params)
        plt = bar(sens_OptStiff_byRange(r),n,yData(n),'FaceColor', 'flat');
        plt.CData = colors(n,:);
    end
    subplotLabel(r) = text(0.025,0.85,charlbl{r},'Units','normalized','FontSize',14);
    xl = get(gca, 'xlim');
    hold on
    plot(gca,xl, [1 1], 'r', 'lineWidth', 1)
    if ((1-ParamRange(r))*100) > 0
        ylabel(strcat({'$+'}, num2str((1-ParamRange(r))*100), {' '}, '\%$'),'Interpreter','latex');
    else
        ylabel(strcat({'$'}, num2str((1-ParamRange(r))*100), {' '}, '\%$'),'Interpreter','latex');
    end
    sens_OptStiff_byRange(r).YLim = [0 6.5];
    grid(sens_OptStiff_byRange(r), 'on')
%     if r ~= 4
        sens_OptStiff_byRange(r).XAxis.Visible = 'off';
%     end
end
% leg = legend(sens_OptStiff_byRange(2),params, 'Location','northeastoutside');
leg = legend(sens_OptStiff_byRange(1),params, 'Location','northoutside','Orientation','horizontal', 'Numcolumns',6);
% xlabel(t, 'Substrate stiffness (pN/nm)', 'Interpreter', 'latex')
ylabel(t, 'Parameter sensitivity', 'Interpreter', 'latex', 'fontsize', 14)
% title(t, 'Optimum stiffness', 'Interpreter', 'latex')

figname = ['Sens_optStiff_allParams_',condition];
nicePlot(f2h, 14, 20, 0.75, savedirectory, figname,saveGraphs, showGraphs)

%% Param sensitivity plot - OptStiff - selectParams

f2h = figure;
t = tiledlayout(4,1);
for r = 1:4
    yData = yData_OptStiff_byRange{r}; % returns a 1x22 matrix 
    yData = yData(:,qlist); % subsetting data for the selected parameters
    sens_OptStiff_byRange(r) = nexttile;
    hold on
    for n = 1:numel(plotParams)
        plt = bar(sens_OptStiff_byRange(r),n,yData(n),'FaceColor', 'flat');
        plt.CData = colors(n,:);
    end
    subplotLabel(r) = text(0.025,0.85,charlbl{r},'Units','normalized','FontSize',14);
    xl = get(gca, 'xlim');
    hold on
    plot(gca,xl, [1 1], 'r', 'lineWidth', 1)
    if ((1-ParamRange(r))*100) > 0
        ylabel(strcat({'$+'}, num2str((1-ParamRange(r))*100), {' '}, '\%$'),'Interpreter','latex');
    else
        ylabel(strcat({'$'}, num2str((1-ParamRange(r))*100), {' '}, '\%$'),'Interpreter','latex');
    end
    sens_OptStiff_byRange(r).YLim = [0 6.5];
    grid(sens_OptStiff_byRange(r), 'on')
%     if r ~= 4
        sens_OptStiff_byRange(r).XAxis.Visible = 'off';
%     end
end
% leg = legend(sens_OptStiff_byRange(2),plotParams, 'Location','northeastoutside');
leg = legend(sens_OptStiff_byRange(1),plotParams, 'Location','northoutside','Orientation','horizontal', 'Numcolumns',6);
% xlabel(t, 'Substrate stiffness (pN/nm)', 'Interpreter', 'latex')
ylabel(t, 'Parameter sensitivity', 'Interpreter', 'latex', 'fontsize', 14)
% title(t, 'Optimum stiffness', 'Interpreter', 'latex')

figname = ['Sens_optStiff_selectParams_',condition];
nicePlot(f2h, 14, 20, 0.75, savedirectory, figname,saveGraphs, showGraphs)

%% Plot sensitivity value on y axis vs stiffness - all params

f1h = figure;

for r = 1%:4
    for p = 1:length(params) %cycle through parameters
        IndexC = strfind(params,params(p));
        q = find(not(cellfun('isempty',IndexC)));
        yData = sensitivity(q).IntMHS(r,:); %  for a given range value and stiffness, extract sensitivities of all parameters of interest.     
        s = subplot(4,6,p);
        plot(k_sub_range, yData)
        set(gca,'xscale','log');
%         set(gca,'xlim',[1,295]);
        title(params(p));
        s.XTick = '';
    end
end

figname = ['sens_vs_stiffness_allParams_',condition];
nicePlot(f1h, 14, 30, 0.6, savedirectory, figname, saveGraphs, showGraphs)
    
%% Plot sensitivity value on y axis vs stiffness - select Params

f1h = figure;

for r = 1%:4
    for p = 1:length(plotParams) %cycle through parameters
        IndexC = strfind(params,plotParams(p));
        q = find(not(cellfun('isempty',IndexC)));
        yData = sensitivity(q).IntMHS(r,:); %  for a given range value and stiffness, extract sensitivities of all parameters of interest.     
        s = subplot(4,2,p);
        plot(k_sub_range, yData)
        set(gca,'xscale','log');
%         set(gca,'xlim',[1,295]);
        title(plotParams(p));
        s.XTick = '';
    end
end
   
figname = ['sens_vs_stiffness_selectParams_',condition];
nicePlot(f1h, 14, 20, 0.9, savedirectory, figname, saveGraphs, showGraphs)

%% Unpacking struct
k_sub_range = [0.01, 0.02:0.02:10, 12.5, 15:5:30, 40, 50, 60, 80, 100, 150, 200, 250, 350, 500, 650, 800, 1000];

vels_unpacked = NaN(numel(params), numel(ParamRange), numel(k_sub_range));
intMHS_unpacked = NaN(numel(params), numel(ParamRange), numel(k_sub_range));

for i = 1:size(sensData,2)
    vels_unpacked(i,:,:) = sensData(i).meanVels;
    intMHS_unpacked(i,:,:) = sensData(i).IntMHS;
end

%% Plot velocity curves - all params 

f1h = figure;
t1 = tiledlayout(6,4);
for p = 1:length(params) %cycle through parameters
    IndexC = strfind(params,params(p));
    q = find(not(cellfun('isempty',IndexC)));
    sens_vels_curves(p)= nexttile;
%     subplotLabel(p) = text(0.4,0.85,charlbl{p},'Units','normalized','FontSize',12);
    hold on
    for j = 1:size(vels_unpacked,2) %cycling through ParamRange values
        plot(k_sub_range, squeeze(vels_unpacked(q,j,:)),'DisplayName',num2str(ParamRange(j)),'LineWidth',1);
    end
    [min1, loc1] = min(squeeze(vels_unpacked(q,1,:)));
    [min2, loc2] = min(squeeze(vels_unpacked(q,4,:)));
    d = [k_sub_range(loc2), min2]-[k_sub_range(loc1), min1];
    quiver(k_sub_range(loc1), min1, d(1), d(2), 0, 'LineWidth', 1.5,'color', 'k', 'HandleVisibility','off');
    set(gca,'xscale','log')
    xlim([0, 1e3])
    xticks([10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3]);
%     xticklabels(['10^-2', '10^-1', '10^0', '10^1', '10^2', '10^3']);
    title(params{q});
end
nexttile(24)
plot(k_sub_range, ones(5,length(k_sub_range))*nan, 'LineWidth', 1.2)
set(gca, 'xtick', [],'visible','off')
leg = legend('Base $- 20\%$','Base $- 10\%$','Base $+ 10\%$', 'Base $+ 20\%$','Base value', 'Location', 'eastoutside');
% leg = legend('Location','northeastoutside');
xlabel(t1, 'Substrate stiffness (pN/nm)','Interpreter','latex', 'fontsize', 14)
ylabel(t1, 'Mean retrograde velocity (nm/s)','Interpreter','latex', 'fontsize', 14)
% title(t1, 'Percentage of integrins in mid and high order species at equlibrium', 'Interpreter', 'latex')

figname = ['Sens_vels_curves_allParams_', condition];
nicePlot(f1h, 14, 30, 0.6, savedirectory, figname ,saveGraphs, showGraphs)

%% Plot velocity curves - select params 

f1h = figure;
t1 = tiledlayout(4,2);
for p = 1:length(plotParams) %cycle through parameters
    IndexC = strfind(params,plotParams(p));
    q = find(not(cellfun('isempty',IndexC)));
    sens_vels_curves(p)= nexttile;
%     subplotLabel(p) = text(0.4,0.85,charlbl{p},'Units','normalized','FontSize',12);
    hold on
    for j = 1:size(vels_unpacked,2) %cycling through ParamRange values
        plot(k_sub_range, squeeze(vels_unpacked(q,j,:)),'DisplayName',num2str(ParamRange(j)),'LineWidth',1);
    end
    [min1, loc1] = min(squeeze(vels_unpacked(q,1,:)));
    [min2, loc2] = min(squeeze(vels_unpacked(q,4,:)));
    d = [k_sub_range(loc2), min2]-[k_sub_range(loc1), min1];
    quiver(k_sub_range(loc1), min1, d(1), d(2), 0, 'LineWidth', 1.5,'color', 'k', 'HandleVisibility','off');
    set(gca,'xscale','log')
    xlim([0, 1e3])
    xticks([10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3]);
%     xticklabels(['10^-2', '10^-1', '10^0', '10^1', '10^2', '10^3']);
    title(params{q});
end
% nexttile(9)
% plot(k_sub_range, ones(5,length(k_sub_range))*nan, 'LineWidth', 1.2)
% set(gca, 'xtick', [],'visible','off')
leg = legend('Base $- 20\%$','Base $- 10\%$','Base $+ 10\%$', 'Base $+ 20\%$','Base value', 'Location', 'eastoutside');
% leg = legend('Location','northeastoutside');
xlabel(t1, 'Substrate stiffness (pN/nm)','Interpreter','latex', 'fontsize', 14)
ylabel(t1, 'Mean retrograde velocity (nm/s)','Interpreter','latex', 'fontsize', 14)
% title(t1, 'Percentage of integrins in mid and high order species at equlibrium', 'Interpreter', 'latex')

figname = ['Sens_vels_curves_selectParams_', condition];
nicePlot(f1h, 14, 19, 0.9, savedirectory, figname,saveGraphs, showGraphs);

%% Plot all Maturation Fraction curves - all params
f2h = figure;
t2 = tiledlayout(6,4);
for p = 1:length(params) %cycle through parameters
    IndexC = strfind(params,params(p));
    q = find(not(cellfun('isempty',IndexC)));
    sens_MF_curves(p)= nexttile;
%     subplotLabel(p) = text(0.4,0.85,charlbl{p},'Units','normalized','FontSize',12);
    hold on
    for j = 1:size(intMHS_unpacked,2)
        plot(k_sub_range, squeeze(intMHS_unpacked(q,j,:)),'DisplayName',num2str(ParamRange(j)),'LineWidth',1);
    end
    set(gca,'xscale','log')
    xlim([0, 1e3])
    ylim([0.05, 0.43])
    xticks([10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3]);
    title(params{q});
end
nexttile(24)
plot(k_sub_range, ones(5,length(k_sub_range))*nan, 'LineWidth', 1.2)
set(gca, 'xtick', [],'visible','off')
leg = legend('Base $- 20\%$','Base $- 10\%$','Base $+ 10\%$', 'Base $+ 20\%$','Base value', 'Location', 'eastoutside');
% leg = legend('Location','northeastoutside');
xlabel(t2, 'Substrate stiffness (pN/nm)','Interpreter','latex', 'fontsize', 14)
ylabel(t2, 'Maturation Fraction','Interpreter','latex', 'fontsize', 14)
% title(t1, 'Percentage of integrins in mid and high order species at equlibrium', 'Interpreter', 'latex')

figname = ['Sens_MF_curves_allParams_', condition];
nicePlot(f2h, 14, 30, 0.6, savedirectory, figname, saveGraphs, showGraphs)

%% Plot all Maturation Fraction curves - select params
f2h = figure;
t2 = tiledlayout(4,2);
for p = 1:length(plotParams) %cycle through parameters
    IndexC = strfind(params,plotParams(p));
    q = find(not(cellfun('isempty',IndexC)));
    sens_MF_curves(p)= nexttile;
%     subplotLabel(p) = text(0.4,0.85,charlbl{p},'Units','normalized','FontSize',12);
    hold on
    for j = 1:size(intMHS_unpacked,2)
        plot(k_sub_range, squeeze(intMHS_unpacked(q,j,:)),'DisplayName',num2str(ParamRange(j)),'LineWidth',1);
    end
    set(gca,'xscale','log')
    xlim([0, 1e3])
    ylim([0.05, 0.43])
    xticks([10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3]);
    title(params{q});
end
% nexttile(24)
% plot(k_sub_range, ones(5,length(k_sub_range))*nan, 'LineWidth', 1.2)
% set(gca, 'xtick', [],'visible','off')
leg = legend('Base $- 20\%$','Base $- 10\%$','Base $+ 10\%$', 'Base $+ 20\%$','Base value', 'Location', 'eastoutside');
% leg = legend('Location','northeastoutside');
xlabel(t2, 'Substrate stiffness (pN/nm)','Interpreter','latex', 'fontsize', 14)
ylabel(t2, 'Maturation Fraction','Interpreter','latex', 'fontsize', 14)
% title(t1, 'Percentage of integrins in mid and high order species at equlibrium', 'Interpreter', 'latex')

figname = ['Sens_MF_curves_selectParams_',condition];
nicePlot(f2h, 14, 19, 0.9, savedirectory, figname,saveGraphs, showGraphs)

%% Load experimental actin velocity and traction force data.

velocity_peyton2005 = readmatrix('C:\Users\kaila\Documents\MSc\Thesis project\Git\Calibration data\cell speed - stiffness - biphasic_ from Peyton_2005 ref126.csv');

velocity_Elosegui2016 = readmatrix('C:\Users\kaila\Documents\MSc\Thesis project\Git\Calibration data\actin flow speed - biphasic_ from Elosegui_artola 2016.csv');
velocity_ChanOdde2008 = readmatrix('C:\Users\kaila\Documents\MSc\Thesis project\Git\Calibration data\Chan&Odde2008_measuredActinVelocity.csv');
traction_Elosegui2016 = readmatrix('C:\Users\kaila\Documents\MSc\Thesis project\Git\Calibration data\cell_traction_biphasic_Elosegui-Artola 2016.csv');

velocity_Elosegui2014_comput = readmatrix('C:\Users\kaila\Documents\MSc\Thesis project\Git\Calibration data\Elosegui-artola 2014 computational velocity.csv');

% converting young's modulus to spring constant
velocity_Elosegui2016(:,1) = (4*pi*(550)/9).*(10^-3).*velocity_Elosegui2016(:,1);
velocity_ChanOdde2008(:,1) = (4*pi*(550)/9).*(10^-3).*velocity_ChanOdde2008(:,1);
velocity_Elosegui2014_comput(:,1) = (4*pi*(550)/9).*(10^-3).*velocity_Elosegui2014_comput(:,1);
velocity_peyton2005(:,1) = (4*pi*(550)/9).*(10^-3).*velocity_peyton2005(:,1);

% converting to nm/s
velocity_peyton2005(:,2) = velocity_peyton2005(:,2)*1000/60; % this is cell migration speed, not actin retrograde speed. But the maximum is around 25 kPa. 

%% Biphasic velocity and force curves
vels_unpacked = NaN(numel(params), numel(ParamRange), numel(k_sub_range));
f_unpacked = NaN(numel(params), numel(ParamRange), numel(k_sub_range));

for i = 1:size(sensData,2)
    vels_unpacked(i,:,:) = sensData(i).meanVels;
    f_unpacked(i,:,:) = sensData(i).meanForces;
end

% % Mean vels data for baseline value of any of the 21 parameters can be
% % used to plot the force vs stiffness biphasic curve. Here I pick the 1st
% % parameter, and the 1st row (ie baseline value row for that param)
meanVels = squeeze(vels_unpacked(1,5,:))';
meanF = squeeze(f_unpacked(1,5,:))';

fh8 = figure;
yyaxis left
hold on
plot(k_sub_range, meanVels, 'LineWidth', 1.2);
% plot(velocity_Elosegui2016(:,1), velocity_Elosegui2016(:,2), '-*', 'DisplayName', 'Elosegui-Artola (2016)','LineWidth', 1.2);
% plot(velocity_ChanOdde2008(:,1), velocity_ChanOdde2008(:,2), '-o', 'DisplayName', 'Chan and Odde (2008)','LineWidth', 1.2);
% plot(velocity_Elosegui2014_comput(:,1), velocity_Elosegui2014_comput(:,2), 'DisplayName', 'Elosegui-Artola (predicted) (2014)','LineWidth', 1.2);
set(gca,'xscale','log')
ylabel('Mean velocity (nm/s)');
hold off
yyaxis right
plot(k_sub_range, meanF, 'LineWidth', 1.2)
% title('Mean actin retrograde velocity and force vs substrate stiffness')
xlabel('Substrate stiffness (pN/nm)');
ylabel('Mean force (pN)');
ylab = get(gca, 'YLabel');
set(ylab,'rotation',-90)
ylab.Position(1) = ylab.Position(1)+ 1100;

figname = ['BiphasicVelsForces',condition];
nicePlot(fh8, 14, 20, 0.65, savedirectory, figname, saveGraphs, showGraphs)

% Biphasic plot with data from previous studies
fh1 = figure;
hold on
plot(k_sub_range, meanVels, 'LineWidth', 1.2, 'DisplayName', 'This study');
% plot(velocity_Elosegui2016(:,1), velocity_Elosegui2016(:,2), '-*', 'DisplayName', 'Elosegui-Artola et al. (2016)','LineWidth', 1.2);
plot(velocity_ChanOdde2008(:,1), velocity_ChanOdde2008(:,2), '-o', 'DisplayName', 'Chan and Odde (2008)','LineWidth', 1.2);
plot(velocity_Elosegui2014_comput(:,1), velocity_Elosegui2014_comput(:,2), 'DisplayName', 'Elosegui-Artola et al. (2014) - predicted','LineWidth', 1.2);
set(gca,'xscale','log')
% leg = legend(sens_IntMHS_byRange(1),plotParams, 'Location','northoutside','Orientation','horizontal');
leg = legend('Location','northoutside','Orientation','horizontal', 'NumColumns', 2);
xlabel('Substrate stiffness (pN/nm)', 'fontsize', 14);
ylabel('Mean velocity (nm/s)', 'fontsize', 14);

figname = ['BiphasicVels_ExptData',condition];
nicePlot(fh1, 14, 20, 0.65, savedirectory, figname, saveGraphs, showGraphs)

% Y = 0.1; %kPa (Young's modulus) range to be tested -> 0.1 to 100 kPa 
% a = 550; %nm (radius of adhesion)
% k_sub = (4*pi*a/9).* Y.*(10^-3); % Converting young's modulus to spring constant k_sub. 10^-3 is multiplied to convert units to pN/nm. 

%% Biphasic maturation fraction curves

MF_unpacked = NaN(numel(params), numel(ParamRange), numel(k_sub_range));

for i = 1:size(sensData,2)
    MF_unpacked(i,:,:) = sensData(i).IntMHS;
end

% % Maturation fraction data for baseline value of any of the 21 parameters can be
% % used to plot the force vs stiffness biphasic curve. Here I pick the 1st
% % parameter, and the 5th row (ie baseline value row for that param)
MatFrac = squeeze(MF_unpacked(1,5,:))';

fh8 = figure;
% yyaxis left
% hold on
plot(k_sub_range, MatFrac, 'LineWidth', 1.2);
set(gca,'xscale','log')
ylabel('Maturation fraction');
xlabel('Substrate stiffness (pN/nm)');

figname = ['MF_',condition];
nicePlot(fh8, 14, 20, 0.65, savedirectory, figname, saveGraphs, showGraphs)

%% Tiled plot of velocity, force and MF vs stiffness

fh1 = figure;
t1 = tiledlayout(1,3);
combPlot{1} = nexttile;
hold on
plot(k_sub_range, meanVels, 'LineWidth', 1.2, 'DisplayName', 'This study');
% plot(velocity_Elosegui2016(:,1), velocity_Elosegui2016(:,2), '-*', 'DisplayName', 'Elosegui-Artola et al. (2016)','LineWidth', 1.2);
plot(velocity_ChanOdde2008(:,1), velocity_ChanOdde2008(:,2), '-o', 'DisplayName', 'Chan and Odde (2008)','LineWidth', 1.2);
plot(velocity_Elosegui2014_comput(:,1), velocity_Elosegui2014_comput(:,2), 'DisplayName', 'Elosegui-Artola et al. (2014) - predicted','LineWidth', 1.2);
set(gca,'xscale','log')
% leg = legend(sens_IntMHS_byRange(1),plotParams, 'Location','northoutside','Orientation','horizontal');
leg = legend('Location','northoutside','Orientation','horizontal', 'NumColumns', 1);
xlabel('Substrate stiffness (pN/nm)', 'fontsize', 11);
ylabel('Mean velocity (nm/s)', 'fontsize', 11);
hold off

combPlot{2} = nexttile;
hold on
yyaxis left
hold on
plot(k_sub_range, meanVels, 'LineWidth', 1.2);
% plot(velocity_Elosegui2016(:,1), velocity_Elosegui2016(:,2), '-*', 'DisplayName', 'Elosegui-Artola (2016)','LineWidth', 1.2);
% plot(velocity_ChanOdde2008(:,1), velocity_ChanOdde2008(:,2), '-o', 'DisplayName', 'Chan and Odde (2008)','LineWidth', 1.2);
% plot(velocity_Elosegui2014_comput(:,1), velocity_Elosegui2014_comput(:,2), 'DisplayName', 'Elosegui-Artola (predicted) (2014)','LineWidth', 1.2);
set(gca,'xscale','log')
ylabel('Mean velocity (nm/s)');
hold off
yyaxis right
plot(k_sub_range, meanF, 'LineWidth', 1.2)
% title('Mean actin retrograde velocity and force vs substrate stiffness')
xlabel('Substrate stiffness (pN/nm)');
ylabel('Mean force (pN)');
ylab = get(gca, 'YLabel');
set(ylab,'rotation',-90)
% ylab.Position(1) = ylab.Position(1)+ 2000;

combPlot{3} = nexttile;
plot(k_sub_range, MatFrac, 'LineWidth', 1.2);
set(gca,'xscale','log')
ylabel('Maturation fraction');
xlabel('Substrate stiffness (pN/nm)');
ylab2 = get(gca,'YLabel');
% ylab2.Position(1) = ylab2.Position(1) + 0.0001;

figname = ['combVelsMF_', condition];
nicePlot(fh1, 10, 19, 0.45, savedirectory, figname, saveGraphs, showGraphs)
