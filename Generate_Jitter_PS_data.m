% Generates the injection delay PS simulation data as shown in Manning et al.
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;

addpath('DCE_Simulation_Functions');

[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;

PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;
vP_range = linspace(SimParam.min_vP,SimParam.max_vP,10)';
PS_fixed = PhysParam.PS_fixed;

N_PS = size(PS_range,1); %range sizes to test
N_vP = size(vP_range,1); %range sizes to test
Delay_ranges = [0 4 8 12]; % Injection delay ranges

%% Sim jitter with Patlak fitting (fast injection, Patlak fit)
% T1 acquisition
acqParam.T1_SNR = 318;
for m = 1:N_PS
    for n = 1:SimParam.N_repetitions
        T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
        T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
        [T1_blood_meas_s(n,m),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
        T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
        T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
        [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
    end
end
% loop through PS and vP values, simulate
for i = 1:size(Delay_ranges,2);
     SimParam.t_start_s = 119 + Delay_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_jitter_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_jitter_fast(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     
     PS_means_jitter_fast(:,i) = mean(PS_fit_jitter_fast,1)'; % mean for each PS
     PS_devs_jitter_fast(:,i) = std(PS_fit_jitter_fast,0,1)'; % standard deviation for each PS
     vP_means_jitter_fast(:,i) = mean(vP_fit_jitter_fast,1)'; % mean for each vP
     vP_devs_jitter_fast(:,i) = std(vP_fit_jitter_fast,0,1)'; % standard deviation for each vP
end
 
%% Sim jitter with Patlak fitting (fast injection, Patlak fit, w/ exclusion fit)
SimParam.SXLfit = 0;
SimParam.NIgnore = max(SimParam.baselineScans) + 3;

for m = 1:N_PS
    for n = 1:SimParam.N_repetitions
        T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
        T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
        [T1_blood_meas_s(n,m),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
        T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
        T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
        [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
    end
end
% loop through PS and vP values, simulate
for i = 1:size(Delay_ranges,2);
     SimParam.t_start_s = 119 + Delay_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_jitter_fast_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_jitter_fast_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     
     PS_means_jitter_fast_exclude(:,i) = mean(PS_fit_jitter_fast_exclude,1)'; % mean for each PS
     PS_devs_jitter_fast_exclude(:,i) = std(PS_fit_jitter_fast_exclude,0,1)'; % standard deviation for each PS
     vP_means_jitter_fast_exclude(:,i) = mean(vP_fit_jitter_fast_exclude,1)'; % mean for each vP
     vP_devs_jitter_fast_exclude(:,i) = std(vP_fit_jitter_fast_exclude,0,1)'; % standard deviation for each vP
end

%% Sim jitter with Patlak fitting (fast injection, NXL fit)
SimParam.SXLfit = 1; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans);

    for m = 1:N_PS
    for n = 1:SimParam.N_repetitions
        T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
        T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
        [T1_blood_meas_s(n,m),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
        T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
        T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
        [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
    end
end
% loop through PS and vP values, simulate
for i = 1:size(Delay_ranges,2);
     SimParam.t_start_s = 119 + Delay_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_jitter_fast_SXL(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_jitter_fast_SXL(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     
     PS_means_jitter_fast_SXL(:,i) = mean(PS_fit_jitter_fast_SXL,1)'; % mean for each PS
     PS_devs_jitter_fast_SXL(:,i) = std(PS_fit_jitter_fast_SXL,0,1)'; % standard deviation for each PS
     vP_means_jitter_fast_SXL(:,i) = mean(vP_fit_jitter_fast_SXL,1)'; % mean for each vP
     vP_devs_jitter_fast_SXL(:,i) = std(vP_fit_jitter_fast_SXL,0,1)'; % standard deviation for each vP
end


%% Sim jitter with Patlak fitting (fast injection, NXL fit, w/ exclusion)
SimParam.SXLfit = 1; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;

for m = 1:N_PS
    for n = 1:SimParam.N_repetitions
        T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
        T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
        [T1_blood_meas_s(n,m),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
        T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
        T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
        [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
    end
end
% loop through PS and vP values, simulate
for i = 1:size(Delay_ranges,2);
    SimParam.t_start_s = 119 + Delay_ranges(i);
    for i_PS = 1:N_PS
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_jitter_fast_SXL_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    for i_vP = 1:N_vP
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
        PhysParam.PS = PS_fixed(1);
        PhysParam.vP = vP_range(i_vP);
        [vP_fit_jitter_fast_SXL_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_jitter_fast_SXL_exclude(:,i) = mean(PS_fit_jitter_fast_SXL_exclude,1)'; % mean for each PS
    PS_devs_jitter_fast_SXL_exclude(:,i) = std(PS_fit_jitter_fast_SXL_exclude,0,1)'; % standard deviation for each PS
    vP_means_jitter_fast_SXL_exclude(:,i) = mean(vP_fit_jitter_fast_SXL_exclude,1)'; % mean for each vP
    vP_devs_jitter_fast_SXL_exclude(:,i) = std(vP_fit_jitter_fast_SXL_exclude,0,1)'; % standard deviation for each vP
end


%% Sim jitter with Patlak fitting (slow injection, Patlak fit)
SimParam.SXLfit = 0;
SimParam.NIgnore = max(SimParam.baselineScans);

SimParam.t_start_s = 0;
SimParam.InjectionRate = 'slow';

 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
 
 % T1 acquisition
 acqParam.T1_SNR = 318;
 for m = 1:N_PS
     for n = 1:SimParam.N_repetitions
         T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
         T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
         [T1_blood_meas_s(n,m),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
         T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
         T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
         [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
     end
 end
 
 % loop through PS and vP values, simulate
 for i = 1:size(Delay_ranges,2);
      SimParam.t_start_s = Delay_ranges(i);
      for i_PS = 1:N_PS
          PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
          PhysParam.vP = vP_fixed(1);
          PhysParam.PS_perMin = PS_range(i_PS);
          [temp, PS_fit_jitter_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      for i_vP = 1:N_vP
          PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
          PhysParam.PS = PS_fixed(1);
          PhysParam.vP = vP_range(i_vP);
          [vP_fit_jitter_slow(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      PS_means_jitter_slow(:,i) = mean(PS_fit_jitter_slow,1)'; % mean for each PS
      PS_devs_jitter_slow(:,i) = std(PS_fit_jitter_slow,0,1)'; % standard deviation for each PS
      vP_means_jitter_slow(:,i) = mean(vP_fit_jitter_slow,1)'; % mean for each vP
      vP_devs_jitter_slow(:,i) = std(vP_fit_jitter_slow,0,1)'; % standard deviation for each vP
 end

  %% Sim jitter with Patlak fitting (slow injection, Patlak fit, w/ exclusion)
      SimParam.SXLfit = 0; % fit enhancements according to SXL method
    SimParam.NIgnore = max(SimParam.baselineScans) + 3;
    
   for m = 1:N_PS
     for n = 1:SimParam.N_repetitions
         T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
         T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
         [T1_blood_meas_s(n,m),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
         T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
         T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
         [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
     end
 end
 
 % loop through PS and vP values, simulate
 for i = 1:size(Delay_ranges,2);
      SimParam.t_start_s = Delay_ranges(i);
      for i_PS = 1:N_PS
          PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
          PhysParam.vP = vP_fixed(1);
          PhysParam.PS_perMin = PS_range(i_PS);
          [temp, PS_fit_jitter_slow_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      for i_vP = 1:N_vP
          PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
          PhysParam.PS = PS_fixed(1);
          PhysParam.vP = vP_range(i_vP);
          [vP_fit_jitter_slow_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      PS_means_jitter_slow_exclude(:,i) = mean(PS_fit_jitter_slow_exclude,1)'; % mean for each PS
      PS_devs_jitter_slow_exclude(:,i) = std(PS_fit_jitter_slow_exclude,0,1)'; % standard deviation for each PS
      vP_means_jitter_slow_exclude(:,i) = mean(vP_fit_jitter_slow_exclude,1)'; % mean for each vP
      vP_devs_jitter_slow_exclude(:,i) = std(vP_fit_jitter_slow_exclude,0,1)'; % standard deviation for each vP
 end
 %% Sim jitter with Patlak fitting (slow injection, NXL fit)
 SimParam.SXLfit = 1; % fit enhancements according to SXL method
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 for m = 1:N_PS
     for n = 1:SimParam.N_repetitions
         T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
         T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
         [T1_blood_meas_s(n,m),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
         T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
         T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
         [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
     end
 end
 
 % loop through PS and vP values, simulate
 for i = 1:size(Delay_ranges,2);
      SimParam.t_start_s = Delay_ranges(i);
      for i_PS = 1:N_PS
          PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
          PhysParam.vP = vP_fixed(1);
          PhysParam.PS_perMin = PS_range(i_PS);
          [temp, PS_fit_jitter_slow_SXL(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      for i_vP = 1:N_vP
          PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
          PhysParam.PS = PS_fixed(1);
          PhysParam.vP = vP_range(i_vP);
          [vP_fit_jitter_slow_SXL(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      PS_means_jitter_slow_SXL(:,i) = mean(PS_fit_jitter_slow_SXL,1)'; % mean for each PS
      PS_devs_jitter_slow_SXL(:,i) = std(PS_fit_jitter_slow_SXL,0,1)'; % standard deviation for each PS
      vP_means_jitter_slow_SXL(:,i) = mean(vP_fit_jitter_slow_SXL,1)'; % mean for each vP
      vP_devs_jitter_slow_SXL(:,i) = std(vP_fit_jitter_slow_SXL,0,1)'; % standard deviation for each vP
 end
 %% Sim jitter with Patlak fitting (slow injection, NXL fit, w/ exclusion)
 SimParam.SXLfit = 1; % fit enhancements according to SXL method
 SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 
 for m = 1:N_PS
     for n = 1:SimParam.N_repetitions
         T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
         T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
         [T1_blood_meas_s(n,m),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
         T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
         T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
         [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
     end
 end
 
 % loop through PS and vP values, simulate
 for i = 1:size(Delay_ranges,2);
      SimParam.t_start_s = Delay_ranges(i);
      for i_PS = 1:N_PS
          PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
          PhysParam.vP = vP_fixed(1);
          PhysParam.PS_perMin = PS_range(i_PS);
          [temp, PS_fit_jitter_slow_SXL_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      for i_vP = 1:N_vP
          PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
          PhysParam.PS = PS_fixed(1);
          PhysParam.vP = vP_range(i_vP);
          [vP_fit_jitter_slow_SXL_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      PS_means_jitter_slow_SXL_exclude(:,i) = mean(PS_fit_jitter_slow_SXL_exclude,1)'; % mean for each PS
      PS_devs_jitter_slow_SXL_exclude(:,i) = std(PS_fit_jitter_slow_SXL_exclude,0,1)'; % standard deviation for each PS
      vP_means_jitter_slow_SXL_exclude(:,i) = mean(vP_fit_jitter_slow_SXL_exclude,1)'; % mean for each vP
      vP_devs_jitter_slow_SXL_exclude(:,i) = std(vP_fit_jitter_slow_SXL_exclude,0,1)'; % standard deviation for each vP
 end
 
%% convert units, save results
PS_range = PS_range * 1e4;
PS_means_jitter_fast = PS_means_jitter_fast * 1e4;
PS_means_jitter_slow = PS_means_jitter_slow * 1e4;
PS_means_jitter_fast_SXL = PS_means_jitter_fast_SXL * 1e4;
PS_means_jitter_slow_SXL = PS_means_jitter_slow_SXL * 1e4;
PS_means_jitter_fast_exclude = PS_means_jitter_fast_exclude * 1e4;
PS_means_jitter_slow_exclude = PS_means_jitter_slow_exclude * 1e4;
PS_means_jitter_fast_SXL_exclude = PS_means_jitter_fast_SXL_exclude * 1e4;
PS_means_jitter_slow_SXL_exclude = PS_means_jitter_slow_SXL_exclude * 1e4;
PS_devs_jitter_fast = PS_devs_jitter_fast * 1e4;
PS_devs_jitter_slow = PS_devs_jitter_slow * 1e4;
PS_devs_jitter_fast_SXL = PS_devs_jitter_fast_SXL * 1e4;
PS_devs_jitter_slow_SXL = PS_devs_jitter_slow_SXL * 1e4;
PS_devs_jitter_fast_exclude = PS_devs_jitter_fast_exclude * 1e4;
PS_devs_jitter_slow_exclude = PS_devs_jitter_slow_exclude * 1e4;
PS_devs_jitter_fast_SXL_exclude = PS_devs_jitter_fast_SXL_exclude * 1e4;
PS_devs_jitter_slow_SXL_exclude = PS_devs_jitter_slow_SXL_exclude * 1e4;

vP_range = vP_range * 1e2;
vP_means_jitter_fast = vP_means_jitter_fast * 1e2;
vP_means_jitter_slow = vP_means_jitter_slow * 1e2;
vP_means_jitter_fast_SXL = vP_means_jitter_fast_SXL * 1e2;
vP_means_jitter_slow_SXL = vP_means_jitter_slow_SXL * 1e2;
vP_means_jitter_fast_exclude = vP_means_jitter_fast_exclude * 1e2;
vP_means_jitter_slow_exclude = vP_means_jitter_slow_exclude * 1e2;
vP_means_jitter_fast_SXL_exclude = vP_means_jitter_fast_SXL_exclude * 1e2;
vP_means_jitter_slow_SXL_exclude = vP_means_jitter_slow_SXL_exclude * 1e2;
vP_devs_jitter_fast = vP_devs_jitter_fast * 1e2;
vP_devs_jitter_slow = vP_devs_jitter_slow * 1e2;
vP_devs_jitter_fast_SXL = vP_devs_jitter_fast_SXL * 1e2;
vP_devs_jitter_slow_SXL = vP_devs_jitter_slow_SXL * 1e2;
vP_devs_jitter_fast_exclude = vP_devs_jitter_fast_exclude * 1e2;
vP_devs_jitter_slow_exclude = vP_devs_jitter_slow_exclude * 1e2;
vP_devs_jitter_fast_SXL_exclude = vP_devs_jitter_fast_SXL_exclude * 1e2;
vP_devs_jitter_slow_SXL_exclude = vP_devs_jitter_slow_SXL_exclude * 1e2;

save('PS_means_jitter','PS_means_jitter_fast','PS_means_jitter_fast_SXL','PS_means_jitter_slow','PS_means_jitter_slow_SXL'...
    ,'PS_means_jitter_fast_exclude','PS_means_jitter_slow_exclude','PS_means_jitter_fast_SXL_exclude','PS_means_jitter_slow_SXL_exclude')
save('PS_devs_jitter','PS_devs_jitter_fast','PS_devs_jitter_fast_SXL','PS_devs_jitter_slow','PS_devs_jitter_slow_SXL'...
    ,'PS_devs_jitter_fast_exclude','PS_devs_jitter_slow_exclude','PS_devs_jitter_fast_SXL_exclude','PS_devs_jitter_slow_SXL_exclude')
save('vP_means_jitter','vP_means_jitter_fast','vP_means_jitter_fast_SXL','vP_means_jitter_slow','vP_means_jitter_slow_SXL'...
    ,'vP_means_jitter_fast_exclude','vP_means_jitter_slow_exclude','vP_means_jitter_fast_SXL_exclude','vP_means_jitter_slow_SXL_exclude')
save('vP_devs_jitter','vP_devs_jitter_fast','vP_devs_jitter_fast_SXL','vP_devs_jitter_slow','vP_devs_jitter_slow_SXL'...
    ,'vP_devs_jitter_fast_exclude','vP_devs_jitter_slow_exclude','vP_devs_jitter_fast_SXL_exclude','vP_devs_jitter_slow_SXL_exclude')


%% Plot figure of results
Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];
Colour4 = [0.4940 0.1840 0.5560 0.5];

figure()

h1=subplot(2,4,1) % fast, FXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_fast(:,1) - PS_range, 1*PS_devs_jitter_fast(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_fast(:,2) - PS_range, 1*PS_devs_jitter_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_fast(:,3) - PS_range, 1*PS_devs_jitter_fast(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_fast(:,4) - PS_range, 1*PS_devs_jitter_fast(:,4),'LineWidth',1.1,'Color',Colour4);
ylabel({'{\bfNo exclusion}'},'FontSize',8);
title('Bolus (FXL fitting)','FontSize',8);
xlim([0 max(PS_range)]);
ylim([-5 5]);
legend({'No delay','+ 4 s','+ 8 s','+ 12 s'},'Location','best','FontSize',8)
legend('boxoff')

subplot(2,4,2) % slow, FXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_slow(:,1) - PS_range, 1*PS_devs_jitter_slow(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_slow(:,2) - PS_range, 1*PS_devs_jitter_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_slow(:,3) - PS_range, 1*PS_devs_jitter_slow(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_slow(:,4) - PS_range, 1*PS_devs_jitter_slow(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([0 max(PS_range)]);
title('Slow inj (FXL fitting)','FontSize',8);
ylim([-1.25 1.6]);

subplot(2,4,3) % fast, NXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_fast_SXL(:,1) - PS_range, 1*PS_devs_jitter_fast_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_fast_SXL(:,2) - PS_range, 1*PS_devs_jitter_fast_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_fast_SXL(:,3) - PS_range, 1*PS_devs_jitter_fast_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_fast_SXL(:,4) - PS_range, 1*PS_devs_jitter_fast_SXL(:,4),'LineWidth',1.1,'Color',Colour4);
title('Bolus (NXL fitting)','FontSize',8);
xlim([0 max(PS_range)]);
ylim([-1.25 1.6]);

subplot(2,4,4) % slow, NXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_slow_SXL(:,1) - PS_range, 1*PS_devs_jitter_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_slow_SXL(:,2) - PS_range, 1*PS_devs_jitter_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_slow_SXL(:,3) - PS_range, 1*PS_devs_jitter_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_slow_SXL(:,4) - PS_range, 1*PS_devs_jitter_slow_SXL(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([0 max(PS_range)]);
title('Slow inj (NXL fitting)','FontSize',8);
ylim([-1.25 1.6]);

h2=subplot(2,4,5) % fast, FXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_fast_exclude(:,1) - PS_range, 1*PS_devs_jitter_fast_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_fast_exclude(:,2) - PS_range, 1*PS_devs_jitter_fast_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_fast_exclude(:,3) - PS_range, 1*PS_devs_jitter_fast_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_fast_exclude(:,4) - PS_range, 1*PS_devs_jitter_fast_exclude(:,4),'LineWidth',1.1,'Color',Colour4);
ylabel({'{\bfw/ exclusion}'},'FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)]);
ylim([-1.25 1.6]);

subplot(2,4,6) % slow, FXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_slow_exclude(:,1) - PS_range, 1*PS_devs_jitter_slow_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_slow_exclude(:,2) - PS_range, 1*PS_devs_jitter_slow_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_slow_exclude(:,3) - PS_range, 1*PS_devs_jitter_slow_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_slow_exclude(:,4) - PS_range, 1*PS_devs_jitter_slow_exclude(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([0 max(PS_range)]);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
ylim([-1.25 1.6]);

subplot(2,4,7) % fast, NXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_fast_SXL_exclude(:,1) - PS_range, 1*PS_devs_jitter_fast_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_fast_SXL_exclude(:,2) - PS_range, 1*PS_devs_jitter_fast_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_fast_SXL_exclude(:,3) - PS_range, 1*PS_devs_jitter_fast_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_fast_SXL_exclude(:,4) - PS_range, 1*PS_devs_jitter_fast_SXL_exclude(:,4),'LineWidth',1.1,'Color',Colour4);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)]);
ylim([-1.25 1.6]);

subplot(2,4,8) % slow, NXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_slow_SXL_exclude(:,1) - PS_range, 1*PS_devs_jitter_slow_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_slow_SXL_exclude(:,2) - PS_range, 1*PS_devs_jitter_slow_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_slow_SXL_exclude(:,3) - PS_range, 1*PS_devs_jitter_slow_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_slow_SXL_exclude(:,4) - PS_range, 1*PS_devs_jitter_slow_SXL_exclude(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([0 max(PS_range)]);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
ylim([-1 1.6]);

p1=get(h1,'position');
p2=get(h2,'position');
height=p1(2)+p1(4)-p2(2);
hx1=axes('position',[0.11 p2(2) p2(3) height],'visible','off');
h_label=ylabel('Fitted {\itPS} error (x10^{-4} min^{-1} )','visible','on');

annotation(figure(1),'textbox',[0.090 0.920 0.05 0.045],'String','(A)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.297 0.920 0.06 0.045],'String','(B)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.503 0.920 0.06 0.045],'String','(C)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.709 0.920 0.06 0.045],'String','(D)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.090 0.445 0.06 0.045],'String','(E)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.297 0.445 0.06 0.045],'String','(F)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.503 0.445 0.06 0.045],'String','(G)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.709 0.445 0.06 0.045],'String','(H)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);

set(gcf, 'units', 'centimeters','PaperPosition', [0 0 20 16]);    % can be bigger than screen
print(gcf, 'Figure_4.png', '-dpng','-r800');
print(gcf, 'Figure_4.tif', '-dtiff','-r800');

%% Figure of vp values (for supplementary material)
figure(2)
h3=subplot(2,4,1) % fast, FXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_fast(:,1) - vP_range, 1*vP_devs_jitter_fast(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_jitter_fast(:,2) - vP_range, 1*vP_devs_jitter_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_jitter_fast(:,3) - vP_range, 1*vP_devs_jitter_fast(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.039, vP_means_jitter_fast(:,4) - vP_range, 1*vP_devs_jitter_fast(:,4),'LineWidth',1.1,'Color',Colour4);
ylabel({'{\bfNo exclusion}'},'FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-1.5 0.5]);
title('Bolus (FXL fitting)','FontSize',8);

subplot(2,4,2) % slow, FXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_slow(:,1) - vP_range, 1*vP_devs_jitter_slow(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_jitter_slow(:,2) - vP_range, 1*vP_devs_jitter_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_jitter_slow(:,3) - vP_range, 1*vP_devs_jitter_slow(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.039, vP_means_jitter_slow(:,4) - vP_range, 1*vP_devs_jitter_slow(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.5]);
title('Slow inj (FXL fitting)','FontSize',8);

subplot(2,4,3) % fast, NXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_fast_SXL(:,1) - vP_range, 1*vP_devs_jitter_fast_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_jitter_fast_SXL(:,2) - vP_range, 1*vP_devs_jitter_fast_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_jitter_fast_SXL(:,3) - vP_range, 1*vP_devs_jitter_fast_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.039, vP_means_jitter_fast_SXL(:,4) - vP_range, 1*vP_devs_jitter_fast_SXL(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.5]);
title('Bolus (NXL fitting)','FontSize',8);

subplot(2,4,4) % slow, NXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_slow_SXL(:,1) - vP_range, 1*vP_devs_jitter_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_jitter_slow_SXL(:,2) - vP_range, 1*vP_devs_jitter_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_jitter_slow_SXL(:,3) - vP_range, 1*vP_devs_jitter_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.039, vP_means_jitter_slow_SXL(:,4) - vP_range, 1*vP_devs_jitter_slow_SXL(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.5]);
title('Slow inj (NXL fitting)','FontSize',8);

h4=subplot(2,4,5) % fast, FXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_fast_exclude(:,1) - vP_range, 1*vP_devs_jitter_fast_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_jitter_fast_exclude(:,2) - vP_range, 1*vP_devs_jitter_fast_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_jitter_fast_exclude(:,3) - vP_range, 1*vP_devs_jitter_fast_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.039, vP_means_jitter_fast_exclude(:,4) - vP_range, 1*vP_devs_jitter_fast_exclude(:,4),'LineWidth',1.1,'Color',Colour4);
ylabel({'{\bfw/ exclusion}'},'FontSize',8);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.5]);
legend({'No delay','+ 4 s','+ 8 s','+ 12 s'},'Location','best','FontSize',8)
legend('boxoff')

subplot(2,4,6) % slow, FXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_slow_exclude(:,1) - vP_range, 1*vP_devs_jitter_slow_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_jitter_slow_exclude(:,2) - vP_range, 1*vP_devs_jitter_slow_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_jitter_slow_exclude(:,3) - vP_range, 1*vP_devs_jitter_slow_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.039, vP_means_jitter_slow_exclude(:,4) - vP_range, 1*vP_devs_jitter_slow_exclude(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([min(vP_range) max(vP_range)+0.026]);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
ylim([-0.5 0.5]);

subplot(2,4,7) % fast, NXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_fast_SXL_exclude(:,1) - vP_range, 1*vP_devs_jitter_fast_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_jitter_fast_SXL_exclude(:,2) - vP_range, 1*vP_devs_jitter_fast_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_jitter_fast_SXL_exclude(:,3) - vP_range, 1*vP_devs_jitter_fast_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.039, vP_means_jitter_fast_SXL_exclude(:,4) - vP_range, 1*vP_devs_jitter_fast_SXL_exclude(:,4),'LineWidth',1.1,'Color',Colour4);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.5]);

subplot(2,4,8) % slow, NXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_slow_SXL_exclude(:,1) - vP_range, 1*vP_devs_jitter_slow_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_jitter_slow_SXL_exclude(:,2) - vP_range, 1*vP_devs_jitter_slow_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_jitter_slow_SXL_exclude(:,3) - vP_range, 1*vP_devs_jitter_slow_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.039, vP_means_jitter_slow_SXL_exclude(:,4) - vP_range, 1*vP_devs_jitter_slow_SXL_exclude(:,4),'LineWidth',1.1,'Color',Colour4);
xlim([min(vP_range) max(vP_range)+0.026]);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
ylim([-0.5 0.5]);



p3=get(h3,'position');
p4=get(h4,'position');
height=p3(2)+p3(4)-p4(2);
hx2=axes('position',[0.11 p4(2) p4(3) height],'visible','off');
h_label=ylabel('Fitted {\itv_p} error (x10^{-2})','visible','on');

set(gcf, 'units', 'centimeters','Position', [5 5 20 16]);

annotation(figure(2),'textbox',[0.090 0.920 0.06 0.045],'String','(A)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(2),'textbox',[0.297 0.920 0.06 0.045],'String','(B)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(2),'textbox',[0.503 0.920 0.06 0.045],'String','(C)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(2),'textbox',[0.709 0.920 0.06 0.045],'String','(D)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(2),'textbox',[0.090 0.445 0.06 0.045],'String','(E)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(2),'textbox',[0.297 0.445 0.06 0.045],'String','(F)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(2),'textbox',[0.503 0.445 0.06 0.045],'String','(G)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(2),'textbox',[0.709 0.445 0.06 0.045],'String','(H)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);

set(gcf, 'units', 'centimeters','PaperPosition', [0 0 20 16]);    % can be bigger than screen
set(gcf,'PaperPositionMode','manual');
print(gcf, 'Supp_Figure_vP_jitter.png', '-dpng','-r800');

