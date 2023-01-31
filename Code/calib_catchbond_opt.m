% Parameter optimization
%% integrin fibronectin catch bond

% Best params for the function A*exp(-b*Force) + C*exp(d*Force) were:
% A = 9.26035213810656;
% b = 0.162850645095768;
% C = 0.000546353153716452;
% d = 0.158421544238933; 

% With a5b1 data
%  Closer param values to other comp papers: 2.51675514872035         0.106855913680857      0.000134394969149903         0.189677012719433
%  or even this :                            2.1909979116314         0.102491397049064      0.000265238563667941         0.170816557905542

% With a5b3 data
% data = readmatrix('C:/Users/kaila/Documents/MSc/Thesis project/Git/Calibration data/Calib_catch-bond_a5b3integrin_ref6.csv');

% a5b1 data from 
data = readmatrix('D:\Thesis project\Master Folder\Calibration data\Calib_catch-bond_a5b1integrin_kong2009.csv');

% Catch-bond 
forceData = data(:,1);
lifetimeData = data(:,2);

func_catch_param_opt = @(x)sseval_catchbond(x,forceData,lifetimeData);


close all
% startingPoints = rand(4,1);
startingPoints = [2+0.1*rand,1+0.1*rand,0.1+0.1*rand,0.1+0.1*rand];
options = optimset('Display','final','MaxIter',10000000,'MaxFunEvals',10000000,'TolFun',1e-8,'TolX',1e-8);%,'PlotFcns',@optimplotfval);
[bestFitParams,fval,exitflag,fitInfo] = fminsearch(func_catch_param_opt, startingPoints);

A = bestFitParams(1);
b = bestFitParams(2);
C = bestFitParams(3);
d = bestFitParams(4);

x = [0:0.01:60];
y1 = 1./(A*exp(-b*x)+C*exp(d*x));
plot(x,y1,'b','DisplayName','Fitted curve');
hold on
plot(forceData,lifetimeData,'*', 'DisplayName','Experimental data')
legend
xlabel('Force (pN)')
ylabel('Lifetime (s)')
bestFitParams
fval
% exitflag
% fitInfo

%% talin unfolding rate as a function of force - talin-rod

% %  best fit params for k_unf = k_unf_UL*exp(k*F) are 
% k_unf_UL = 1.53690747659748;
% k =  0.0501825624410996;

% using this, k_unf for 12 pN was 2.8, very close to what they got from
% this curve in the paper with this data - https://www.science.org/doi/full/10.1126/science.1162912
data = readmatrix('C:/Users/kaila/Documents/MSc/Thesis project/Git/Calibration data/Calib_slip-bond_talin-unfold_ref11_Rio2009.csv');

% talin unfolding
forceData = data(:,1);
rData = data(:,2);

func_slip_param_opt = @(x)sseval_slip(x,forceData,rData);


close all
% startingPoints = rand(4,1);
startingPoints = [0.01*rand,0.01*rand];
options = optimset('Display','final','MaxIter',10000000,'MaxFunEvals',10000000,'TolFun',1e-8,'TolX',1e-8);%,'PlotFcns',@optimplotfval);
[bestFitParams,fval,exitflag,fitInfo] = fminsearch(func_slip_param_opt, startingPoints);

k_unf_UL = bestFitParams(1);
k = bestFitParams(2);

% k_unf_UL = 0.018;
% k = bestFitParams(2)+0.21;

x = [0:0.01:50];
y2 = k_unf_UL*exp(k*x);
plot(x,y2,'b','DisplayName','Fitted curve');
hold on
plot(forceData,rData,'*', 'DisplayName','Experimental data')
bestFitParams
fval
% exitflag
% fitInfo

y3 = func_tal_unf(x,2,maxRate);
plot(x,y3,'r','DisplayName','Original parameters');
legend
xlabel('Force (pN)')
ylabel('Rate (s^-1)')

%% talin R3 unfolding from Tapia-Rojo et al 2020 
% In model 11 the fitted parameters from del rio et al 2009 was used. 
% k_unf_UL = 1.53690747659748;
% k =  0.0501825624410996;

% Talin_unf_tapiaRojo2020
%  k_unf_UL = 0.00592416109336837
%  k = 1.58006965404131;
data = readmatrix('C:/Users/kaila/Documents/MSc/Thesis project/Git/Calibration data/Talin_unf_tapiaRojo2020.csv');

% talin unfolding
forceData = data(:,1);
rData = data(:,2);

func_slip_param_opt = @(x)sseval_slip(x,forceData,rData);


close all
% startingPoints = rand(4,1);
startingPoints = [0.1*rand,0.1*rand];
options = optimset('Display','final','MaxIter',10000000,'MaxFunEvals',10000000,'TolFun',1e-8,'TolX',1e-8);%,'PlotFcns',@optimplotfval);
[bestFitParams,fval,exitflag,fitInfo] = fminsearch(func_slip_param_opt, startingPoints);

k_unf_UL = bestFitParams(1);
k = bestFitParams(2);

x = [0:0.01:20];
% y2 = k_unf_UL*exp(k*x);
y2 = log(k_unf_UL)+(k*x);
y3 = exp(y2);
plot(x,y2,'b','DisplayName','Fitted curve');
hold on
plot(forceData,rData,'*', 'DisplayName','Experimental data')
bestFitParams
fval
% exitflag
% fitInfo

%% Talin-unfolding yao et al (mechanical response of talin)
% Talin_unf_Yao2016 (mechanical response of talin)

%  k_unf_UL = 0.0174375902587974 
%  k = 1.38515568984156;
data = readmatrix('C:/Users/kaila/Documents/MSc/Thesis project/Git/Calibration data/Talin_unf_Yao2016 (mechanical response of talin).csv');

% talin unfolding
forceData = data(:,1);
rData = data(:,2);

func_slip_param_opt = @(x)sseval_slip(x,forceData,rData);


close all
% startingPoints = rand(4,1);
startingPoints = [0.1*rand,0.1*rand];
options = optimset('Display','final','MaxIter',10000000,'MaxFunEvals',10000000,'TolFun',1e-8,'TolX',1e-8);%,'PlotFcns',@optimplotfval);
[bestFitParams,fval,exitflag,fitInfo] = fminsearch(func_slip_param_opt, startingPoints);

k_unf_UL = bestFitParams(1);
k = bestFitParams(2);

x = [0:0.01:20];
% y2 = k_unf_UL*exp(k*x);
y2 = log(k_unf_UL)+(k*x);
y3 = exp(y2);
plot(x,y2,'b','DisplayName','Fitted curve');
hold on
plot(forceData,rData,'*', 'DisplayName','Experimental data')
bestFitParams
fval
% exitflag
% fitInfo
tal_unf = @(k_unf_UL,k,F) k_unf_UL*exp(k*F);

%% Percentage of maturing and non maturing adhesions
% 38.3 percent of NAs that form mature into focal complexes. (Paper - 100% of G2
% mature into focal complexes, of which 32% mature into focal adhesions) - https://elifesciences.org/articles/66151


data = readmatrix('C:/Users/kaila/Documents/MSc/Thesis project/Git/Calibration data/Number of NAs formed and matured_han et al.csv');
data(:,1) = [1;2;3;4;5];
bar(data(:,2))
hold on
plot(data(:,1),data(:,2),'*');

totalNAs = sum(data(:,2));
percMaturing = data(2,2)/totalNAs;
percNonMaturing = sum(data([1,3,4,5],2))/totalNAs;

%% NA assembly disassembly data
NA_assemblyData = readmatrix('C:/Users/kaila/Documents/MSc/Thesis project/Git/Calibration data/NA_assembly_disassembly_choi et al.csv');
NA_assemblyData_norm = NA_assemblyData;
NA_assemblyData_norm(:,2) = normalize(NA_assemblyData_norm(:,2),'range');

%%
% func_slip_param_opt = @(x)sseval_slip(x,tData,data);
% 
% close all
% % startingPoints = rand(4,1);
% startingPoints = [rand];
% options = optimset('Display','final','MaxIter',10000000,'MaxFunEvals',10000000,'TolFun',1e-8,'TolX',1e-8);%,'PlotFcns',@optimplotfval);
% [bestFitParams,fval,exitflag,fitInfo] = fminsearch(func_slip_param_opt, startingPoints);
% 
% k_unf_UL = bestFitParams(1);
% x = [0:0.01:600];
% % y2 = k_unf_UL*exp(k*x);
% y2 = exp(-k_unf_UL*x);
% y3 = exp(y2);
% plot(x,y2,'b','DisplayName','Fitted curve');
% hold on
% plot(tData,data,'DisplayName','Experimental data')
% legend
% 
% bestFitParams
fval
