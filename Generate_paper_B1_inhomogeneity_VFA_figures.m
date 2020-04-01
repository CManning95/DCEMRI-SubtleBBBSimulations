% Generates the PS and vP simulation figures for B1 inhomogeneity effects on T1 acquisition using variable flip angle as shown in Manning et al.
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
acqParam.T1_acq_method = 'Accurate';  % T1 acquisition method ('VFA' or 'HIFI')
acqParam.isFit = [1 1 1]; % which acqusitions to fit (fit all for HIFI, last 3 for VFA)
acqParam.TR_s = [0.0054 0.0054 0.0054]; % repetition times for T1 acqusition
acqParam.FA_nom_rads = [2 5 12] *2*(pi/360); % nominal flip angles for T1 acquisition
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
acqParam.isIR = [0 0 0]; % indicates which are IR-SPGR
acqParam.TI_s = [NaN NaN NaN]; % Inversion times for T1 acquisition (for HIFI)
acqParam.PECentre = [NaN NaN NaN]; % indicates time of centre of k-space
acqParam.NReadout = [160 160 160]; % number of readout pulses (Siemens - number of slices)
acqParam.NTry = 1; % fitting attempts

%% Generate variable water exchange figures

%derive additional parameters
DCESeqParam.NPoints = round(DCESeqParam.t_acq_s/DCESeqParam.t_res_sample_s); % number of sample points

% ranges of PS and vP to test
PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;

PS_fixed = PhysParam.PS_fixed;
vP_range = linspace(SimParam.min_vP,SimParam.max_vP,10)'+1e-8;

%range sizes to test
N_PS = size(PS_range,1);
N_vP = size(vP_range,1);

 % Generate variable temporal jitter PS and vP - fast injection (no exclude)
 SimParam.InjectionRate = 'fast';
  SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
  SimParam.NIgnore = max(SimParam.baselineScans);

% Simulate T1 acquisiton
acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1; % Accurate DCE
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_Accurate_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,1) = mean(PS_fit_Accurate_fast,1)'; % mean for each PS for k = 1
PS_devs_VFA_fast(:,1) = std(PS_fit_Accurate_fast,0,1)'; % standard deviation at k = 1

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_Accurate_fast(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_fast(:,1) = mean(vP_fit_Accurate_fast,1)'; % mean for each vP for k = 1
vP_devs_VFA_fast(:,1) = std(vP_fit_Accurate_fast,0,1)'; % standard deviation for k = 1

% Simulate T1 acquisiton
acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Inaccurate DCE
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_correctedT1_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,2) = mean(PS_fit_correctedT1_fast,1)'; % add mean for each PS for k = 1.2
PS_devs_VFA_fast(:,2) = std(PS_fit_correctedT1_fast,0,1)'; % add standard deviation for k = 1.2

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_correctedT1_fast(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_fast(:,2) = mean(vP_fit_correctedT1_fast,1)'; % add mean for each vP for k = 1.2
vP_devs_VFA_fast(:,2) = std(vP_fit_correctedT1_fast,0,1)'; % add standard deviation for k = 1.2

% Simulate T1 acquisiton
acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Inaccurate DCE
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_NoCorrection_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,3) = mean(PS_fit_NoCorrection_fast,1)'; % add mean for each PS for k = 0.8
PS_devs_VFA_fast(:,3) = std(PS_fit_NoCorrection_fast,0,1)'; % add standard deviation for k = 0.8

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_NoCorrection_fast(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_fast(:,3) = mean(vP_fit_NoCorrection_fast,1)'; % add mean for each vP for k = 0.8
vP_devs_VFA_fast(:,3) = std(vP_fit_NoCorrection_fast,0,1)'; % add standard deviation for k = 0.8

 % Generate variable flow PS and vP - fast injection (exclude)
  SimParam.InjectionRate = 'fast';
   SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans) + 3;
  
% Simulate T1 acquisiton
acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1; % Accurate DCE
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

  for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_accurate_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_exclude(:,1) = mean(PS_fit_accurate_exclude,1)'; % mean for each PS for k = 1
PS_devs_VFA_exclude(:,1) = std(PS_fit_accurate_exclude,0,1)'; % standard deviation for k = 1

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_accurate_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_exclude(:,1) = mean(vP_fit_accurate_exclude,1)'; % mean for each vP for k = 1
vP_devs_VFA_exclude(:,1) = std(vP_fit_accurate_exclude,0,1)'; % standard deviation for k = 1

DCESeqParam.FA_error = 1.2;
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
% Simulate T1 acquisiton
acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Accurate DCE
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_correctedT1_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_exclude(:,2) = mean(PS_fit_correctedT1_exclude,1)'; % mean for each PS for k = 1.2
PS_devs_VFA_exclude(:,2) = std(PS_fit_correctedT1_exclude,0,1)'; % standard deviation for k = 1.2

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_correctedT1_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_exclude(:,2) = mean(vP_fit_correctedT1_exclude,1)'; % add mean for each vP for k = 1.2
vP_devs_VFA_exclude(:,2) = std(vP_fit_correctedT1_exclude,0,1)'; % add standard deviation for k = 1.2

% Simulate T1 acquisiton
acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Accurate DCE
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_NoCorrection_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_exclude(:,3) = mean(PS_fit_NoCorrection_exclude,1)'; % add mean for each PS for k = 0.8
PS_devs_VFA_exclude(:,3) = std(PS_fit_NoCorrection_exclude,0,1)'; % add standard deviation for k = 0.8

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_NoCorrection_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_exclude(:,3) = mean(vP_fit_NoCorrection_exclude,1)'; % add mean for each vP for k = 0.8
vP_devs_VFA_exclude(:,3) = std(vP_fit_NoCorrection_exclude,0,1)'; % add standard deviation for k = 0.8

 % Generate variable flow PS and vP - slow injection
 SimParam.InjectionRate = 'slow';
  SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 69; % number of time points

% Simulate T1 acquisiton
acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1; % Accurate DCE
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_accurate_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_slow(:,1) = mean(PS_fit_accurate_slow,1)'; % mean for each PS for k = 1
PS_devs_VFA_slow(:,1) = std(PS_fit_accurate_slow,0,1)'; % standard deviation for k = 1

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_accurate_slow(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_slow(:,1) = mean(vP_fit_accurate_slow,1)'; % mean for each vP for k = 1
vP_devs_VFA_slow(:,1) = std(vP_fit_accurate_slow,0,1)'; % standard deviation for k = 1

% Simulate T1 acquisiton
acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Accurate DCE
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_correctedT1_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_slow(:,2) = mean(PS_fit_correctedT1_slow,1)'; % add mean for each PS for k = 1.2
PS_devs_VFA_slow(:,2) = std(PS_fit_correctedT1_slow,0,1)'; % add standard deviation for k = 1.2

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_correctedT1_slow(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_slow(:,2) = mean(vP_fit_correctedT1_slow,1)'; % add mean for each vP for k = 1.2
vP_devs_VFA_slow(:,2) = std(vP_fit_correctedT1_slow,0,1)'; % add standard deviation for k = 1.2

% Simulate T1 acquisiton
acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Accurate DCE
acqParam.FA_true_rads = DCESeqParam.FA_error * acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,acqParam,acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,acqParam,acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_NoCorrection_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_slow(:,3) = mean(PS_fit_NoCorrection_slow,1)'; % add mean for each PS for k = 0.8
PS_devs_VFA_slow(:,3) = std(PS_fit_NoCorrection_slow,0,1)'; % add standard deviation for k = 0.8

for i_vP = 1:N_vP
    PhysParam.PS_perMin = PS_fixed(1);
    PhysParam.vP = vP_range(i_vP);
    [vP_fit_NoCorrection_slow(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end

vP_means_VFA_slow(:,3) = mean(vP_fit_NoCorrection_slow,1)'; % add mean for each vP for k = 0.8
vP_devs_VFA_slow(:,3) = std(vP_fit_NoCorrection_slow,0,1)'; % add standard deviation for k = 0.8

% Change scale of PS and vP values (for graphing)
PS_range = PS_range * 1e4;
PS_means_VFA_fast = PS_means_VFA_fast * 1e4;
PS_means_VFA_exclude = PS_means_VFA_exclude * 1e4;
PS_means_VFA_slow = PS_means_VFA_slow * 1e4;
PS_devs_VFA_fast = PS_devs_VFA_fast * 1e4;
PS_devs_VFA_exclude = PS_devs_VFA_exclude * 1e4;
PS_devs_VFA_slow = PS_devs_VFA_slow * 1e4;

vP_range = vP_range * 1e3;
vP_means_VFA_fast = vP_means_VFA_fast * 1e3;
vP_means_VFA_exclude = vP_means_VFA_exclude * 1e3;
vP_means_VFA_slow = vP_means_VFA_slow * 1e3;
vP_devs_VFA_fast = vP_devs_VFA_fast * 1e3;
vP_devs_VFA_exclude = vP_devs_VFA_exclude * 1e3;
vP_devs_VFA_slow = vP_devs_VFA_slow * 1e3;

% Generate graphs

%Set colours and opacity of plots
Colour1  = [0 0.447 0.741];
Colour2 = [0.85 0.325 0.098];
Colour3 = [0.929 0.694 0.125];

figure(2)
subplot(2,3,1)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_fast(:,1) - PS_range, 1*PS_devs_VFA_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_fast(:,2) - PS_range, 1*PS_devs_VFA_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_fast(:,3) - PS_range, 1*PS_devs_VFA_fast(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True {\it PS} (x10^{-4} min^{-1} )');
ylabel('fitted {\it PS} error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-3.2 3.2]);
legend({'Accurate T1 and DCE','Corrected T1','Uncorrected T1'},'Location','southwest')
legend('boxoff')

ax = gca;
ax.FontSize = 9;

subplot(2,3,2)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_exclude(:,1) - PS_range, 1*PS_devs_VFA_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_exclude(:,2) - PS_range, 1*PS_devs_VFA_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_exclude(:,3) - PS_range, 1*PS_devs_VFA_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True {\it PS} (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 9;

subplot(2,3,3)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_slow(:,1) - PS_range, 1*PS_devs_VFA_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_slow(:,2) - PS_range, 1*PS_devs_VFA_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_slow(:,3) - PS_range, 1*PS_devs_VFA_slow(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True {\it PS} (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 9;

subplot(2,3,4)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_VFA_fast(:,1) - vP_range, 1*vP_devs_VFA_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.015, vP_means_VFA_fast(:,2) - vP_range, 1*vP_devs_VFA_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.03, vP_means_VFA_fast(:,3) - vP_range, 1*vP_devs_VFA_fast(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True {\it v_p} (x10^{-3})');
ylabel('fitted {\it v_p} error (x10^{-3})');
xlim([0 max(vP_range)]);
ylim([-6.2 6.2]);

ax = gca;
ax.FontSize = 9;

subplot(2,3,5)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_VFA_exclude(:,1) - vP_range, 1*vP_devs_VFA_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.015, vP_means_VFA_exclude(:,2) - vP_range, 1*vP_devs_VFA_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.03, vP_means_VFA_exclude(:,3) - vP_range, 1*vP_devs_VFA_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True {\it v_p} (x10^{-3})');
xlim([0 max(vP_range)]);
ylim([-3.1 3.1]);

ax = gca;
ax.FontSize = 9;

subplot(2,3,6)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_VFA_slow(:,1) - vP_range, 1*vP_devs_VFA_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.015, vP_means_VFA_slow(:,2) - vP_range, 1*vP_devs_VFA_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.03, vP_means_VFA_slow(:,3) - vP_range, 1*vP_devs_VFA_slow(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True {\it v_p} (x10^{-3})'); 
xlim([0 max(vP_range)]);
ylim([-3.1 3.1]);

ax = gca;
ax.FontSize = 9;

set(gcf,'units','centimeters','position',[5,5,21.0,14.85]);

% Annotate figure markers
annotation(figure(2),'textbox',[0.1 0.927 0.05 0.045],'String',{'a)'},'LineStyle','none','FitBoxToText','off','fontweight','bold');
annotation(figure(2),'textbox',[0.379 0.927 0.06 0.045],'String','b)','LineStyle','none','FitBoxToText','off','fontweight','bold');
annotation(figure(2),'textbox',[0.656 0.927 0.06 0.045],'String','c)','LineStyle','none','FitBoxToText','off','fontweight','bold');
annotation(figure(2),'textbox',[0.1 0.453 0.05 0.045],'String','d)','LineStyle','none','FitBoxToText','off','fontweight','bold');
annotation(figure(2),'textbox',[0.379 0.453 0.06 0.045],'String','e)','LineStyle','none','FitBoxToText','off','fontweight','bold');
annotation(figure(2),'textbox',[0.656 0.453 0.06 0.045],'String','f)','LineStyle','none','FitBoxToText','off','fontweight','bold');