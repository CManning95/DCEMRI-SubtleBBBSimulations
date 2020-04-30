% Generates the PS simulation figures as shown in Manning et al.
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
SimParam.venous_delay_s = 6; % delay of measured VIF compared to AIF
SimParam.t_start_s = 198; % injection delay in seconds
SimParam.InjectionRate = 'slow'; % Speed of injection (slow or fast)
SimParam.syn_model = '2CXM'; % model to simulate contrast leakage
SimParam.water_exch_model = 'FXL'; % water exchange model to generate signals
SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal

% Acquisition parameters (for determining T1)
% - Variable Flip Angle (VFA), or HIFI (accurate)
T1acqParam.T1_acq_method = 'Accurate';  
T1acqParam.isFit = [1 1 1]; % which acquisitions to fit
T1acqParam.TR_s = [0.0054 0.0054 0.0054]; % repetition times for T1 acqusition
T1acqParam.FA_nom_rads = [2 5 12] *2*(pi/360); % nominal flip angles for T1 acquisition
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
T1acqParam.isIR = [0 0 0]; % indicates which are IR-SPGR
T1acqParam.TI_s = [NaN NaN NaN]; % Inversion times for T1 acquisition (for HIFI)
T1acqParam.PECentre = [NaN NaN NaN]; % indicates time of centre of k-space
T1acqParam.NReadout = [160 160 160]; % number of readout pulses (Siemens - number of slices)
T1acqParam.NTry = 1; % fitting attempts

%% Generate variable flow/injection delay/water exchange figures

% Simulate T1 acquisiton
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);

DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg; % measured flip angle is same as nominal

%derive additional parameters
DCESeqParam.NPoints = round(DCESeqParam.t_acq_s/DCESeqParam.t_res_sample_s); % number of sample points
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

% ranges of PS and vP to test
PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;

%range sizes to test
N_PS = size(PS_range,1);

% Variable ranges to test
Fp_ranges = [18 9 4.5]; % Plasma flow ranges
Delay_ranges = [0 6 12]; % Injection delay ranges
kbe_ranges = [2.75 5.5 11.0]; % kbe ranges
T1_diff_ranges = [1 0.8 1.2];

 % Generate variable Fp PS graphs
 SimParam.InjectionRate = 'fast';
 SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 for i = 1:size(Fp_ranges,2) % Generate variable flow data
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_fast(:,i) = mean(PS_fit_Fp_fast,1)'; % mean for each PS - flow
     PS_devs_Fp_fast(:,i) = std(PS_fit_Fp_fast,0,1)'; % standard deviation for each PS - flow
 end
 
 SimParam.InjectionRate = 'fast';
 SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans) + 3;

 for i = 1:size(Fp_ranges,2) % Generate variable flow data
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_exclude(:,i) = mean(PS_fit_Fp_exclude,1)'; % mean for each PS - flow
     PS_devs_Fp_exclude(:,i) = std(PS_fit_Fp_exclude,0,1)'; % standard deviation for each PS - flow
 end 
 SimParam.InjectionRate = 'slow';
 SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 69; % number of time points
 
 for i = 1:size(Fp_ranges,2)
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_slow(:,i) = mean(PS_fit_Fp_slow,1)'; % mean for each PS at high flow
     PS_devs_Fp_slow(:,i) = std(PS_fit_Fp_slow,0,1)'; % standard deviation for each PS at high flow
 end
  
 PhysParam.FP_mlPer100gPerMin = 18; % reset flow rate
 % Generate jitter PS
  SimParam.InjectionRate = 'fast';
  SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
  SimParam.NIgnore = max(SimParam.baselineScans);
  
for i = 1:size(Delay_ranges,2);
     SimParam.t_start_s = 198 + Delay_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_jitter_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

  PS_means_jitter_fast(:,i) = mean(PS_fit_jitter_fast,1)'; % mean for each PS - jitter
  PS_devs_jitter_fast(:,i) = std(PS_fit_jitter_fast,0,1)'; % standard deviation for each PS - jitter
end

  SimParam.InjectionRate = 'fast';
  SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
  SimParam.NIgnore = max(SimParam.baselineScans) + 3;
  
 for i = 1:size(Delay_ranges,2);
     SimParam.t_start_s = 198 + Delay_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_jitter_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

  PS_means_jitter_exclude(:,i) = mean(PS_fit_jitter_exclude,1)'; % mean for each PS - jitter
  PS_devs_jitter_exclude(:,i) = std(PS_fit_jitter_exclude,0,1)'; % standard deviation for each PS - jitter
 end
  
 SimParam.InjectionRate = 'slow';
 SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 69; % number of time points
 
   for i = 1:size(Delay_ranges,2);
      SimParam.t_start_s = Delay_ranges(i);
      for i_PS = 1:N_PS
          PhysParam.vP = vP_fixed(1);
          PhysParam.PS_perMin = PS_range(i_PS);
          [temp, PS_fit_jitter_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      
      PS_means_jitter_slow(:,i) = mean(PS_fit_jitter_slow,1)'; % mean for each PS at 198s delay
      PS_devs_jitter_slow(:,i) = std(PS_fit_jitter_slow,0,1)'; % standard deviation for each PS at 198s delay
   end
  
   SimParam.t_start_s = 198; % reset injection delay
   
% variable water exchange effects PS
  SimParam.InjectionRate = 'fast';
  SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
  SimParam.NIgnore = max(SimParam.baselineScans);
  
  SimParam.water_exch_model = 'FXL';
 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_FXL_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_H2O_fast(:,1) = mean(PS_fit_FXL_fast,1)'; % mean for each PS for FXL
PS_devs_H2O_fast(:,1) = std(PS_fit_FXL_fast,0,1)'; % standard deviation for FXL
 
SimParam.water_exch_model = '2S1XA';
for i = 1:size(kbe_ranges,2);
    PhysParam.kbe_perS = kbe_ranges(i);
    for i_PS = 1:N_PS
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_2S1X_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_H2O_fast(:,i+1) = mean(PS_fit_2S1X_fast,1)'; % add mean for each PS for 2S1X
    PS_devs_H2O_fast(:,i+1) = std(PS_fit_2S1X_fast,0,1)'; % add standard deviation for 2S1X
end
  SimParam.InjectionRate = 'fast';
  SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
  SimParam.NIgnore = max(SimParam.baselineScans) + 3;
  
 SimParam.water_exch_model = 'FXL';
 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_FXL_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_H2O_exclude(:,1) = mean(PS_fit_FXL_exclude,1)'; % mean for each PS for FXL
PS_devs_H2O_exclude(:,1) = std(PS_fit_FXL_exclude,0,1)'; % standard deviation for FXL
 
SimParam.water_exch_model = '2S1XA';
for i = 1:size(kbe_ranges,2);
    PhysParam.kbe_perS = kbe_ranges(i);
    for i_PS = 1:N_PS
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_2S1X_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_H2O_exclude(:,i+1) = mean(PS_fit_2S1X_exclude,1)'; % add mean for each PS for 2S1X
    PS_devs_H2O_exclude(:,i+1) = std(PS_fit_2S1X_exclude,0,1)'; % add standard deviation for 2S1X
end
  
SimParam.InjectionRate = 'slow';
 SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 69; % number of time points
 
 SimParam.water_exch_model = 'FXL';
  for i_PS = 1:N_PS
      PhysParam.vP = vP_fixed(1);
      PhysParam.PS_perMin = PS_range(i_PS);
      [temp, PS_fit_FXL_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
  PS_means_H2O_slow(:,1) = mean(PS_fit_FXL_slow,1)'; % mean for each PS for FXL
  PS_devs_H2O_slow(:,1) = std(PS_fit_FXL_slow,0,1)'; % standard deviation for FXL
  
  SimParam.water_exch_model = '2S1XA';
  for i = 1:size(kbe_ranges,2);
      PhysParam.kbe_perS = kbe_ranges(i);
      for i_PS = 1:N_PS
          PhysParam.vP = vP_fixed(1);
          PhysParam.PS_perMin = PS_range(i_PS);
          [temp, PS_fit_2S1X_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      
      PS_means_H2O_slow(:,i+1) = mean(PS_fit_2S1X_slow,1)'; % add mean for each PS for 2S1X
      PS_devs_H2O_slow(:,i+1) = std(PS_fit_2S1X_slow,0,1)'; % add standard deviation for 2S1X
      
  end

  SimParam.water_exch_model = 'FXL'; % reset water exchange model
  
%% Generate B1 inhomogeneity/assumed T1 figures
% B1 inhomogeneity figures
% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1; % Accurate DCE
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

% B2 inhomogeneity figures (fast injection, no exclude)
SimParam.InjectionRate = 'fast';
SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
SimParam.NIgnore = max(SimParam.baselineScans);

 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_Accurate_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,1) = mean(PS_fit_Accurate_fast,1)'; % mean for each PS for k = 1
PS_devs_VFA_fast(:,1) = std(PS_fit_Accurate_fast,0,1)'; % standard deviation at k = 1

% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Inaccurate DCE
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_correctedT1_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,2) = mean(PS_fit_correctedT1_fast,1)'; % add mean for each PS for k = 1.2
PS_devs_VFA_fast(:,2) = std(PS_fit_correctedT1_fast,0,1)'; % add standard deviation for k = 1.2

% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Inaccurate DCE
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_NoCorrection_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,3) = mean(PS_fit_NoCorrection_fast,1)'; % add mean for each PS for k = 0.8
PS_devs_VFA_fast(:,3) = std(PS_fit_NoCorrection_fast,0,1)'; % add standard deviation for k = 0.8

 % Generate variable flow PS and vP - fast injection (exclude)
  SimParam.InjectionRate = 'fast';
   SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans) + 3;
  
% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1; % Accurate DCE
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

  for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_accurate_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_exclude(:,1) = mean(PS_fit_accurate_exclude,1)'; % mean for each PS for k = 1
PS_devs_VFA_exclude(:,1) = std(PS_fit_accurate_exclude,0,1)'; % standard deviation for k = 1

DCESeqParam.FA_error = 1.2;
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Accurate DCE
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_correctedT1_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_exclude(:,2) = mean(PS_fit_correctedT1_exclude,1)'; % mean for each PS for k = 1.2
PS_devs_VFA_exclude(:,2) = std(PS_fit_correctedT1_exclude,0,1)'; % standard deviation for k = 1.2

% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Accurate DCE
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS 
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_NoCorrection_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_exclude(:,3) = mean(PS_fit_NoCorrection_exclude,1)'; % add mean for each PS for k = 0.8
PS_devs_VFA_exclude(:,3) = std(PS_fit_NoCorrection_exclude,0,1)'; % add standard deviation for k = 0.8

 % Generate B1 inhomogeneity figures - slow injection
 SimParam.InjectionRate = 'slow';
  SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 69; % number of time points

% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1; % Accurate DCE
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_accurate_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_slow(:,1) = mean(PS_fit_accurate_slow,1)'; % mean for each PS for k = 1
PS_devs_VFA_slow(:,1) = std(PS_fit_accurate_slow,0,1)'; % standard deviation for k = 1

% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'Accurate'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Accurate DCE
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_correctedT1_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_slow(:,2) = mean(PS_fit_correctedT1_slow,1)'; % add mean for each PS for k = 1.2
PS_devs_VFA_slow(:,2) = std(PS_fit_correctedT1_slow,0,1)'; % add standard deviation for k = 1.2

% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.2; % Accurate DCE
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_NoCorrection_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_slow(:,3) = mean(PS_fit_NoCorrection_slow,1)'; % add mean for each PS for k = 0.8
PS_devs_VFA_slow(:,3) = std(PS_fit_NoCorrection_slow,0,1)'; % add standard deviation for k = 0.8

 % Assumed T1 figures
 % Generate PS figure - fast injection (no exclude)
 SimParam.InjectionRate = 'fast';
 SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 for i = 1:size(T1_diff_ranges,2)
     PhysParam.T1_blood_meas_s = PhysParam.T10_blood_s * T1_diff_ranges(i);
     PhysParam.T1_tissue_meas_s = PhysParam.T10_tissue_s * T1_diff_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_T1assumed_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     PS_means_T1assumed_fast(:,i) = mean(PS_fit_T1assumed_fast,1)'; % mean for each PS
     PS_devs_T1assumed_fast(:,i) = std(PS_fit_T1assumed_fast,0,1)'; % standard deviation for each PS     
 end
 % Generate PS figure - fast injection (exclude)
  SimParam.InjectionRate = 'fast';
  SimParam.baselineScans = [3:5]; % datapoints to use for calculating base signal
  SimParam.NIgnore = max(SimParam.baselineScans) + 3;
  
  for i = 1:size(T1_diff_ranges,2)
     PhysParam.T1_blood_meas_s = PhysParam.T10_blood_s * T1_diff_ranges(i);
     PhysParam.T1_tissue_meas_s = PhysParam.T10_tissue_s * T1_diff_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_T1assumed_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     PS_means_T1assumed_exclude(:,i) = mean(PS_fit_T1assumed_exclude,1)'; % mean for each PS
     PS_devs_T1assumed_exclude(:,i) = std(PS_fit_T1assumed_exclude,0,1)'; % standard deviation for each PS
     
  end
 
 % Generate PS figure - slow injection
 SimParam.InjectionRate = 'slow';
 SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 18.49; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 69; % number of time points

  for i = 1:size(T1_diff_ranges,2)
     PhysParam.T1_blood_meas_s = PhysParam.T10_blood_s * T1_diff_ranges(i);
     PhysParam.T1_tissue_meas_s = PhysParam.T10_tissue_s * T1_diff_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_T1assumed_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     PS_means_T1assumed_slow(:,i) = mean(PS_fit_T1assumed_slow,1)'; % mean for each PS at high flow
     PS_devs_T1assumed_slow(:,i) = std(PS_fit_T1assumed_slow,0,1)'; % standard deviation for each PS at high flow
     
  end

%% Generate graphs

% Change scale of PS values (for graphing)
PS_range = PS_range * 1e4;
PS_means_Fp_fast = PS_means_Fp_fast * 1e4;
PS_means_Fp_exclude = PS_means_Fp_exclude * 1e4;
PS_means_Fp_slow = PS_means_Fp_slow * 1e4;
PS_devs_Fp_fast = PS_devs_Fp_fast * 1e4;
PS_devs_Fp_exclude = PS_devs_Fp_exclude * 1e4;
PS_devs_Fp_slow = PS_devs_Fp_slow * 1e4;

PS_means_jitter_fast = PS_means_jitter_fast * 1e4;
PS_means_jitter_exclude = PS_means_jitter_exclude * 1e4;
PS_means_jitter_slow = PS_means_jitter_slow * 1e4;
PS_devs_jitter_fast = PS_devs_jitter_fast * 1e4;
PS_devs_jitter_exclude = PS_devs_jitter_exclude * 1e4;
PS_devs_jitter_slow = PS_devs_jitter_slow * 1e4;

PS_means_H2O_fast = PS_means_H2O_fast * 1e4;
PS_means_H2O_exclude = PS_means_H2O_exclude * 1e4;
PS_means_H2O_slow = PS_means_H2O_slow * 1e4;
PS_devs_H2O_fast = PS_devs_H2O_fast * 1e4;
PS_devs_H2O_exclude = PS_devs_H2O_exclude * 1e4;
PS_devs_H2O_slow = PS_devs_H2O_slow * 1e4;

PS_means_VFA_fast = PS_means_VFA_fast * 1e4;
PS_means_VFA_exclude = PS_means_VFA_exclude * 1e4;
PS_means_VFA_slow = PS_means_VFA_slow * 1e4;
PS_devs_VFA_fast = PS_devs_VFA_fast * 1e4;
PS_devs_VFA_exclude = PS_devs_VFA_exclude * 1e4;
PS_devs_VFA_slow = PS_devs_VFA_slow * 1e4;

PS_means_T1assumed_fast = PS_means_T1assumed_fast * 1e4;
PS_means_T1assumed_exclude = PS_means_T1assumed_exclude * 1e4;
PS_means_T1assumed_slow = PS_means_T1assumed_slow * 1e4;
PS_devs_T1assumed_fast = PS_devs_T1assumed_fast * 1e4;
PS_devs_T1assumed_exclude = PS_devs_T1assumed_exclude * 1e4;
PS_devs_T1assumed_slow = PS_devs_T1assumed_slow * 1e4;
