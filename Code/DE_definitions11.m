function [DE] = DE_definitions11(concs, var_rates, consts)

% % Reactions
% 
% % Precomplex formation
% 
% R1 =	(k1f)*int*tal*vinc - (k1r)*Pcomp;
% 
% % Seed formation
% R2 = 	(k2f)*Pcomp^25 - (k2r)*S1;
% 
% % Cluster formation
% R3 =	(k3f)*S1^2 - (k3r)*C1;
% 
% % Seed activation 
% R4 = 	(k4f)*S1 - (k4r)*S1a ;
% R5 =	(k5f)*S2 - (k5r)*S2a;
% R6 =	(k6f)*S3 - (k6r)*S3a;
% 
% % Seed reinforcement
% R7 =	(k7f)*S1a*vinc^25 - (k7r)*S2a;
% R8 =	(k8f)*S2a*vinc^25 - (k8r)*S3a;
% 
% % Cluster activation
% R9 =	(k9f)*C1 - (k9r)*C1a;
% R10 = (k10f)*C2 - (k10r)*C2a;
% R11 =	(k11f)*C3 - (k11r)*C3a;
% 
% % Cluster reinforcement
% R12 = (k12f)*C1a*vinc^50 - (k12r)*C2a;
% R13 =	(k13f)*C2a*vinc^50 - (k13r)*C3a;
% 
% % Cluster formation from active seed
% R14 = (k14f)*S1a^2 - (k14r)*C1a;
% R15 =	(k15f)*S2a^2 - (k15r)*C2a;
% R16 = (k16f)*S3a^2 - (k16r)*C3a;
% 
% % Talin refolding - Inactive Seed
% R17 =	(k17f)*S3 - (k17r)*S2*vinc^25;
% R18 =	(k18f)*S2 - (k18r)*S1*vinc^25;
% 
% % Talin refolding - Inactive cluster
% R19 =	(k19f)*C3 - (k19r)*C2*vinc^50;
% R20 =	(k20f)*C2 - (k20r)*C1*vinc^50;
% 
% % Cluster breakdown to seed
% R21 =	(k21f)*C3 - (k21r)*S3^2;
% R22 =	(k22f)*C2 - (k22r)*S2^2;

% % Signal molecule 

% R23: Sig --> inSig, MM kinetics. d[inSig]/dt = V_max * [Sig]/(K_m + [Sig])
% d[Sig]/dt = V_max * [Sig]/(K_m + [Sig])


% % Concentrations
int = concs(1);
tal = concs(2);
vinc = concs(3);
Pcomp = concs(4);
S1 = concs(5);
S2 = concs(6);
S3 = concs(7);
S1a = concs(8);
S2a = concs(9);
S3a = concs(10);
C1 = concs(11);
C2 = concs(12);
C3 = concs(13);
C1a = concs(14);
C2a = concs(15);
C3a = concs(16);
Sig = concs(17);
inSig = concs(18);
% Sig2 = concs(19);

% Time/Force dependent rate constants
k4r = var_rates(1);
k5r = var_rates(2);
k6r = var_rates(3);
k7f = var_rates(4);
k8f = var_rates(5);
k9r = var_rates(6);
k10r = var_rates(7);
k11r = var_rates(8);
k12f = var_rates(9);
k13f = var_rates(10);

%  Fixed rate constants
k1f = consts(1);
k1r = consts(2);
k2f = consts(3);
k2r = consts(4);
k3f = consts(5);
k3r = consts(6);
k4f = consts(7);
k5f = consts(8);
k6f = consts(9);
k7r = consts(10);
k8r = consts(11);
k9f = consts(12);
k10f = consts(13);
k11f = consts(14);
k12r = consts(15);
k13r = consts(16);
k14f = consts(17);
k14r = consts(18);
k15f = consts(19);
k15r = consts(20);
k16f = consts(21);
k16r = consts(22);
k17f = consts(23); %var (1/f)?
k17r = consts(24); %0
k18f = consts(25); %var (1/f)?
k18r = consts(26); %0
k19f = consts(27); %var (1/f)?
k19r = consts(28); %0
k20f = consts(29); %var (1/f)?
k20r = consts(30); %0
k21f = consts(31); %var (1/f)? 
k21r = consts(32); %0
k22f = consts(33); %var (1/f)? 
k22r = consts(34); %0
k23_vmax = consts(35); 
k23_km = consts(36);
k23f = 0.0397;
%%
%  Fixed rate constants
% k1f = 0;%consts(1);
% k1r = 0;%consts(2);
% k2f = 0;%consts(3);
% k2r = 0;%consts(4);
% k3f = 0;%consts(5);
% k3r = 0;%consts(6);
% k4f = 0;%consts(7);
% k4r = 0;
% k5f = 0;%consts(8);
% k5r = 0;
% k6f = 0;%consts(9);
% k6r = 0;
% k7f = 0;
% k7r = 0;%consts(10);
% k8f = 0;
% k8r = 0;%consts(11);
% k9f = 0;%consts(12);
% k9r = 0;
% k10f = 0;%consts(13);
% k10r = 0;
% k11f = 0;%consts(14);
% k11r = 0;
% k12f = 0;
% k12r = 0;%consts(15);
% k13f = 0;
% k13r = 0;%consts(16);
% k14f = 0;%consts(17);
% k14r = 0;%consts(18);
% k15f = 0;%consts(19);
% k15r = 0;%consts(20);
% k16f = 0;%consts(21);
% k16r = 0;%consts(22);
% k17f = 0;%consts(23); %var (1/f)?
% k17r = 0;%consts(24); %0
% k18f = 0;%consts(25); %var (1/f)?
% k18r = 0;%consts(26); %0
% k19f = 0;%consts(27); %var (1/f)?
% k19r = 0;%consts(28); %0
% k20f = 0;%consts(29); %var (1/f)?
% k20r = 0;%consts(30); %0
% k21f = 0;%consts(31); %var (1/f)? 
% k21r = 0;%consts(32); %0
% k22f = 0;%consts(33); %var (1/f)? 
% k22r = 0;%consts(34); %0
% k23_vmax = 0;%consts(35); 
% k23_km = 0;%consts(36);

%%

% Reactions
order1 = 2;
order2 = 2;

% Precomplex formation

R1 =	(k1f)*int*tal*vinc - (k1r)*Pcomp;

% Seed formation
R2 = 	(k2f)*Pcomp^order1 - (k2r)*S1;

% Cluster formation
R3 =	(k3f)*S1^2 - (k3r)*C1;

% Seed activation 
R4 = 	(k4f)*S1 - (k4r)*S1a ;
R5 =	(k5f)*S2 - (k5r)*S2a;
R6 =	(k6f)*S3 - (k6r)*S3a;

% Seed reinforcement
R7 =	(k7f)*S1a*vinc^order2 - (k7r)*S2a; %order1
R8 =	(k8f)*S2a*vinc^order2 - (k8r)*S3a;

% Cluster activation
R9 =	(k9f)*C1 - (k9r)*C1a;
R10 =   (k10f)*C2 - (k10r)*C2a;
R11 =	(k11f)*C3 - (k11r)*C3a;

% Cluster reinforcement
R12 =   (k12f)*C1a*vinc^order2 - (k12r)*C2a;
R13 =	(k13f)*C2a*vinc^order2 - (k13r)*C3a;

% Cluster formation from active seed
R14 =   (k14f)*S1a^2 - (k14r)*C1a;
R15 =	(k15f)*S2a^2 - (k15r)*C2a;
R16 =   (k16f)*S3a^2 - (k16r)*C3a;

% Talin refolding - Inactive Seed
R17 =	(k17f)*S3 - (k17r)*S2*vinc^order2; %order1
R18 =	(k18f)*S2 - (k18r)*S1*vinc^order2;

% Talin refolding - Inactive cluster
R19 =	(k19f)*C3 - (k19r)*C2*vinc^order2;
R20 =	(k20f)*C2 - (k20r)*C1*vinc^order2;

% Inactive cluster breakdown to seed
R21 =	(k21f)*C3 - (k21r)*S3^2;
R22 =	(k22f)*C2 - (k22r)*S2^2;

% Signal molecule 
R23 = k23_vmax * (Sig/(k23_km + Sig)); 
R24 = k23f*Sig;

% Differential equations
dInt_dt = -R1;

dTal_dt = -R1;

dVinc_dt = -R1 ...
           -25*R7 ...
           -25*R8 ...
           -50*R12 ...
           -50*R13 ...
           +25*R17 ...
           +25*R18 ...
           +50*R19 ...
           +50*R20;

dPcomp_dt = +R1 ...
            -25*R2;
       
dS1_dt = +R2 ...
         -2*R3 ...
         -R4 ...
         +R18;
     
dS2_dt = -R5 ...
         +R17 ...
         -R18 ...
         +2*R22;
     
dS3_dt = -R6 ...
         -R17 ...
         +2*R21;
     
dS1a_dt = +R4 ...
          -R7 ...
          -2*R14;
      
dS2a_dt = +R5 ...
          +R7 ...
          -R8 ...
          -2*R15;
      
dS3a_dt = +R6 ...
          +R8 ...
          -2*R16;       
       
dC1_dt = +R3 ...
         -R9 ...
         +R20;
     
dC2_dt = -R10 ...
         +R19 ...
         -R20 ...
         -R22;
     
dC3_dt = -R11 ...
         -R19 ...
         -R21;
     
dC1a_dt = +R9 ...
          -R12 ...
          +R14;
      
dC2a_dt = +R10 ...
          +R12 ...
          -R13 ...
          +R15;
      
dC3a_dt = +R11 ...
          +R13 ...
          +R16;
      
dSig_dt = -R23;

dinSig_dt = +R23; 
      
% dSig2 = -R24;

DE =  [dInt_dt;
       dTal_dt;
       dVinc_dt;
       dPcomp_dt;
       dS1_dt;
       dS2_dt;
       dS3_dt;
       dS1a_dt;
       dS2a_dt;
       dS3a_dt;
       dC1_dt;
       dC2_dt;
       dC3_dt;
       dC1a_dt;
       dC2a_dt;
       dC3a_dt;
       dSig_dt;
       dinSig_dt];
%        dSig2];
end 


