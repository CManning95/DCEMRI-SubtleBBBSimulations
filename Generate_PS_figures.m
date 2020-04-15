% Generates the PS and vP simulation figures as shown in Manning et al.
% (2020) Slow injection paper
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;

addpath('DCE_Simulation_Functions');

% Select default parameters

%Physiological parameters
PhysParam.Hct = 0.42; % Haematocrit
PhysParam.vE = 0.2; % EES volume fraction
PhysParam.FP_mlPer100gPerMin = 18; % Plasma flow
PhysParam.T10_blood_s = 1.901; % baseline T1 blood
PhysParam.T10_tissue_s = 0.917; % baseline T1 tissue (NAWM)
PhysParam.T2s0_blood_s = 0.191; % baseline T2* blood
PhysParam.T2s0_tissue_s = 0.050; % baseline T2* tissue
PhysParam.S0_blood = 100; % baseline signal blood
PhysParam.S0_tissue = 100; % baseline signal tissue
PhysParam.kbe_perS = 2.5; % blood - EES water exchange rate
PhysParam.kie_perS = 1.7; % intracellular space - EES exchange rate
PhysParam.vP_fixed = 0.0058; % plasma volume fraction
PhysParam.PS_fixed = 2.96 * 1e-4; % surface area - permeability product

% Sequence parameters
DCESeqParam.t_acq_s = 1268; % acquisition duration
DCESeqParam.t_res_sample_s = 39.62; % sample temporal resolution
DCESeqParam.TR_s = 1e-3*3.4; % repetition time
DCESeqParam.TE_s = 1e-3*1.7; % echo time
DCESeqParam.r1_per_mM_per_s = 5.0; % T1 relaxivity of contrast agent
DCESeqParam.r2_per_mM_per_s = 7.1; % T2 relaxivity of contrast agent
DCESeqParam.FA_nom_deg = 15; % nominal flip angle
DCESeqParam.FA_error = 1; % k value (flip angle error)

% Simulation parameters
SimParam.N_repetitions = 1000; % repetitions at each PS or vP (to quantify effects of noise)
SimParam.t_res_full_s = 0.1; % generated temporal resolution of simulations
SimParam.NIgnore = 0; % number of post-contrast points to exclude (always uses 3 points for baseline signal)
SimParam.SNR = 164; % signal to noise ratio
SimParam.drift_pctPerMin = 0; % signal drift
SimParam.min_PS = 0 * 1e-4; % minimum PS to test vP
SimParam.max_PS = 5 * 1e-4; % maximum PS to test vP
SimParam.min_vP = 0; % minimum vP to test PS
SimParam.max_vP = 0.01; % maximum vP to test PS
SimParam.venous_delay_s = 6; % delay of measured VIF compared to AIF
SimParam.t_start_s = 198; % injection delay in seconds
SimParam.InjectionRate = 'slow'; % Speed of injection (slow or fast)
SimParam.syn_model = '2CXM'; % model to simulate contrast leakage
SimParam.water_exch_model = 'FXL'; % water exchange model to generate signals
SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal

% Acquisition parameters (for determining T1)
% - Variable Flip Angle (VFA), or HIFI (accurate)
acqParam.T1_acq_method = 'Accurate';  % T1 acquisition method