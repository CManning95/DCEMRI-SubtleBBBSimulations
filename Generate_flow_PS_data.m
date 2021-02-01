% Generates the plasma flow PS simulation data as shown in Manning et al.
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;
addpath('DCE_Simulation_Functions');
    
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;


%% Generate variable flow/injection delay sims
% ranges of PS and vP to test
PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;
vP_range = linspace(SimParam.min_vP,SimParam.max_vP,10)';
PS_fixed = PhysParam.PS_fixed;

N_PS = size(PS_range,1); %range sizes to test
N_vP = size(vP_range,1); %range sizes to test
Fp_ranges = [11 8.25 5.5]; % Plasma flow ranges

%% Sim flows with Patlak fitting (fast injection, Patlak fit)
 SimParam.InjectionRate = 'fast';
 SimParam.SXLfit = 0;
 
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
 for i = 1:size(Fp_ranges,2) % Generate variable flow data
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_Fp_fast(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_fast(:,i) = mean(PS_fit_Fp_fast,1)'; % mean for each PS
     PS_devs_Fp_fast(:,i) = std(PS_fit_Fp_fast,0,1)'; % standard deviation for each PS
     
     vP_means_Fp_fast(:,i) = mean(vP_fit_Fp_fast,1)'; % mean for each vP
     vP_devs_Fp_fast(:,i) = std(vP_fit_Fp_fast,0,1)'; % standard deviation for each vP
 end
 
 %% Sim flows with Patlak fitting (fast injection, Patlak fit, w/ exclusion)
  SimParam.InjectionRate = 'fast';
 SimParam.SXLfit = 0;
 SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 
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
 for i = 1:size(Fp_ranges,2) % Generate variable flow data
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_fast_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_Fp_fast_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_fast_exclude(:,i) = mean(PS_fit_Fp_fast_exclude,1)'; % mean for each PS
     PS_devs_Fp_fast_exclude(:,i) = std(PS_fit_Fp_fast_exclude,0,1)'; % standard deviation for each PS
     
     vP_means_Fp_fast_exclude(:,i) = mean(vP_fit_Fp_fast_exclude,1)'; % mean for each vP
     vP_devs_Fp_fast_exclude(:,i) = std(vP_fit_Fp_fast_exclude,0,1)'; % standard deviation for each vP
 end
 
 %% Sim flows with Patlak fitting (fast injection, SXL fit)
 SimParam.InjectionRate = 'fast';
  SimParam.SXLfit = 1;
 SimParam.NIgnore = max(SimParam.baselineScans);
  
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
 for i = 1:size(Fp_ranges,2) % Generate variable flow data
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_fast_SXL(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_Fp_fast_SXL(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
% loop through PS and vP values, simulate
     PS_means_Fp_fast_SXL(:,i) = mean(PS_fit_Fp_fast_SXL,1)'; % mean for each PS
     PS_devs_Fp_fast_SXL(:,i) = std(PS_fit_Fp_fast_SXL,0,1)'; % standard deviation for each PS
     
     vP_means_Fp_fast_SXL(:,i) = mean(vP_fit_Fp_fast_SXL,1)'; % mean for each vP
     vP_devs_Fp_fast_SXL(:,i) = std(vP_fit_Fp_fast_SXL,0,1)'; % standard deviation for each vP
 end 
 
  %% Sim flows with Patlak fitting (fast injection, SXL fit, w/ exclusion)
  SimParam.InjectionRate = 'fast';
  SimParam.SXLfit = 1;
 SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 
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
 for i = 1:size(Fp_ranges,2) % Generate variable flow data
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_SXL_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_Fp_SXL_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
% loop through PS and vP values, simulate
     PS_means_Fp_fast_SXL_exclude(:,i) = mean(PS_fit_Fp_SXL_exclude,1)'; % mean for each PS
     PS_devs_Fp_fast_SXL_exclude(:,i) = std(PS_fit_Fp_SXL_exclude,0,1)'; % standard deviation for each PS
     
     vP_means_Fp_fast_SXL_exclude(:,i) = mean(vP_fit_Fp_SXL_exclude,1)'; % mean for each vP
     vP_devs_Fp_fast_SXL_exclude(:,i) = std(vP_fit_Fp_SXL_exclude,0,1)'; % standard deviation for each vP
 end 
 
 %% Sim flows with Patlak fitting (slow injection, Patlak fit)
 SimParam.InjectionRate = 'slow';
 SimParam.t_start_s = 0;
 SimParam.SXLfit = 0;
 SimParam.NIgnore = max(SimParam.baselineScans);
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
for i = 1:size(Fp_ranges,2)
    PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
    for i_PS = 1:N_PS
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_Fp_slow(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_slow(:,i) = mean(PS_fit_Fp_slow,1)'; % mean for each PS
     PS_devs_Fp_slow(:,i) = std(PS_fit_Fp_slow,0,1)'; % standard deviation for each PS
     
     vP_means_Fp_slow(:,i) = mean(vP_fit_Fp_slow,1)'; % mean for each vP
     vP_devs_Fp_slow(:,i) = std(vP_fit_Fp_slow,0,1)'; % standard deviation for each vP
 end
 
%% Sim flows with Patlak fitting (slow injection, Patlak fit, w/ exclusion)
 SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 SimParam.SXLfit = 0;

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
 for i = 1:size(Fp_ranges,2)
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_slow_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_Fp_slow_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_slow_exclude(:,i) = mean(PS_fit_Fp_slow_exclude,1)'; % mean for each PS
     PS_devs_Fp_slow_exclude(:,i) = std(PS_fit_Fp_slow_exclude,0,1)'; % standard deviation for each PS 
     
     vP_means_Fp_slow_exclude(:,i) = mean(vP_fit_Fp_slow_exclude,1)'; % mean for each vP
     vP_devs_Fp_slow_exclude(:,i) = std(vP_fit_Fp_slow_exclude,0,1)'; % standard deviation for each vP
 end
 
 %% Sim flows with Patlak fitting (slow injection, SXL fit)
  SimParam.NIgnore = max(SimParam.baselineScans);
  SimParam.SXLfit = 1;
  
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
for i = 1:size(Fp_ranges,2)
    PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
    for i_PS = 1:N_PS
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_slow_SXL(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_Fp_slow_SXL(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_slow_SXL(:,i) = mean(PS_fit_Fp_slow_SXL,1)'; % mean for each PS
     PS_devs_Fp_slow_SXL(:,i) = std(PS_fit_Fp_slow_SXL,0,1)'; % standard deviation for each PS
     
     vP_means_Fp_slow_SXL(:,i) = mean(vP_fit_Fp_slow_SXL,1)'; % mean for each vP
     vP_devs_Fp_slow_SXL(:,i) = std(vP_fit_Fp_slow_SXL,0,1)'; % standard deviation for each vP
end
 
 %% Sim flows with Patlak fitting (slow injection, SXL fit, w/ exclusion)
  SimParam.NIgnore = max(SimParam.baselineScans) + 3;
  SimParam.SXLfit = 1;
   
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
 for i = 1:size(Fp_ranges,2)
     PhysParam.FP_mlPer100gPerMin = Fp_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_Fp_slow_SXL_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_Fp_slow_SXL_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end

     PS_means_Fp_slow_SXL_exclude(:,i) = mean(PS_fit_Fp_slow_SXL_exclude,1)'; % mean for each PS
     PS_devs_Fp_slow_SXL_exclude(:,i) = std(PS_fit_Fp_slow_SXL_exclude,0,1)'; % standard deviation for each PS 
     
     vP_means_Fp_slow_SXL_exclude(:,i) = mean(vP_fit_Fp_slow_SXL_exclude,1)'; % mean for each vP
     vP_devs_Fp_slow_SXL_exclude(:,i) = std(vP_fit_Fp_slow_SXL_exclude,0,1)'; % standard deviation for each vP
 end
 
%% Generate graphs

% Change scale of PS values (for graphing)
PS_range = PS_range * 1e4;
PS_means_Fp_fast = PS_means_Fp_fast * 1e4;
PS_means_Fp_fast_exclude = PS_means_Fp_fast_exclude * 1e4;
PS_means_Fp_fast_SXL = PS_means_Fp_fast_SXL * 1e4;
PS_means_Fp_fast_SXL_exclude = PS_means_Fp_fast_SXL_exclude * 1e4;
PS_means_Fp_slow = PS_means_Fp_slow * 1e4;
PS_means_Fp_slow_exclude = PS_means_Fp_slow_exclude * 1e4;
PS_means_Fp_slow_SXL = PS_means_Fp_slow_SXL * 1e4;
PS_means_Fp_slow_SXL_exclude = PS_means_Fp_slow_SXL_exclude * 1e4;
PS_devs_Fp_fast = PS_devs_Fp_fast * 1e4;
PS_devs_Fp_fast_exclude = PS_devs_Fp_fast_exclude * 1e4;
PS_devs_Fp_slow = PS_devs_Fp_slow * 1e4;
PS_devs_Fp_slow_exclude = PS_devs_Fp_slow_exclude * 1e4;
PS_devs_Fp_fast_SXL = PS_devs_Fp_fast_SXL * 1e4;
PS_devs_Fp_fast_SXL_exclude = PS_devs_Fp_fast_SXL_exclude * 1e4;
PS_devs_Fp_slow_SXL = PS_devs_Fp_slow_SXL * 1e4;
PS_devs_Fp_slow_SXL_exclude = PS_devs_Fp_slow_SXL_exclude * 1e4;

vP_range = vP_range * 1e2;
vP_means_Fp_fast = vP_means_Fp_fast * 1e2;
vP_means_Fp_fast_exclude = vP_means_Fp_fast_exclude * 1e2;
vP_means_Fp_slow = vP_means_Fp_slow * 1e2;
vP_means_Fp_slow_exclude = vP_means_Fp_slow_exclude * 1e2;
vP_means_Fp_fast_SXL = vP_means_Fp_fast_SXL * 1e2;
vP_means_Fp_fast_SXL_exclude = vP_means_Fp_fast_SXL_exclude * 1e2;
vP_means_Fp_slow_SXL = vP_means_Fp_slow_SXL * 1e2;
vP_means_Fp_slow_SXL_exclude = vP_means_Fp_slow_SXL_exclude * 1e2;
vP_devs_Fp_fast = vP_devs_Fp_fast * 1e2;
vP_devs_Fp_fast_exclude = vP_devs_Fp_fast_exclude * 1e2;
vP_devs_Fp_slow = vP_devs_Fp_slow * 1e2;
vP_devs_Fp_slow_exclude = vP_devs_Fp_slow_exclude * 1e2;
vP_devs_Fp_fast_SXL = vP_devs_Fp_fast_SXL * 1e2;
vP_devs_Fp_fast_SXL_exclude = vP_devs_Fp_fast_SXL_exclude * 1e2;
vP_devs_Fp_slow_SXL = vP_devs_Fp_slow_SXL * 1e2;
vP_devs_Fp_slow_SXL_exclude = vP_devs_Fp_slow_SXL_exclude * 1e2;


save('PS_means_Fp','PS_means_Fp_fast','PS_means_Fp_fast_exclude','PS_means_Fp_slow','PS_means_Fp_slow_exclude'...
    ,'PS_means_Fp_fast_SXL','PS_means_Fp_fast_SXL_exclude','PS_means_Fp_slow_SXL','PS_means_Fp_slow_SXL_exclude')
save('PS_devs_Fp','PS_devs_Fp_fast','PS_devs_Fp_fast_exclude','PS_devs_Fp_slow','PS_devs_Fp_slow_exclude'...
    ,'PS_devs_Fp_fast_SXL','PS_devs_Fp_fast_SXL_exclude','PS_devs_Fp_slow_SXL','PS_devs_Fp_slow_SXL_exclude')
save('vP_means_Fp','vP_means_Fp_fast','vP_means_Fp_fast_exclude','vP_means_Fp_slow','vP_means_Fp_slow_exclude'...
    ,'vP_means_Fp_fast_SXL','vP_means_Fp_fast_SXL_exclude','vP_means_Fp_slow_SXL','vP_means_Fp_slow_SXL_exclude')
save('vP_devs_Fp','vP_devs_Fp_fast','vP_devs_Fp_fast_exclude','vP_devs_Fp_slow','vP_devs_Fp_slow_exclude'...
    ,'vP_devs_Fp_fast_SXL','vP_devs_Fp_fast_SXL_exclude','vP_devs_Fp_slow_SXL','vP_devs_Fp_slow_SXL_exclude')


%%
Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];

figure()
h1=subplot(4,4,1) % fast, FXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_fast(:,1) - PS_range, 1*PS_devs_Fp_fast(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_fast(:,2) - PS_range, 1*PS_devs_Fp_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_Fp_fast(:,3) - PS_range, 1*PS_devs_Fp_fast(:,3),'LineWidth',1.1,'Color',Colour3);
ylabel({'{\bfNo exclusion}'},'FontSize',8);
%xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
title('Bolus (FXL fitting)','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-0.8 5]);
legend({'{\itF_p} = 11 ml 100g^{-1}min^{-1}','{\itF_p} = 8.25 ml 100g^{-1}min^{-1}','{\itF_p} = 5.5 ml 100g^{-1}min^{-1}'},'Position',[0.172 0.815 0.081 0.0207],'FontSize',4.9)
legend('boxoff')

subplot(4,4,2) % slow, FXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_slow(:,1) - PS_range, 1*PS_devs_Fp_slow(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_slow(:,2) - PS_range, 1*PS_devs_Fp_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_Fp_slow(:,3) - PS_range, 1*PS_devs_Fp_slow(:,3),'LineWidth',1.1,'Color',Colour3);
title('Slow (FXL fitting)','FontSize',8);
%xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-1.6 1.6]);

subplot(4,4,3) % fast, NXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_fast_SXL(:,1) - PS_range, 1*PS_devs_Fp_fast_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_fast_SXL(:,2) - PS_range, 1*PS_devs_Fp_fast_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_Fp_fast_SXL(:,3) - PS_range, 1*PS_devs_Fp_fast_SXL(:,3),'LineWidth',1.1,'Color',Colour3);
%xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
title('Bolus (NXL fitting)','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-1.6 1.6]);

subplot(4,4,4) % slow, NXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_slow_SXL(:,1) - PS_range, 1*PS_devs_Fp_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_slow_SXL(:,2) - PS_range, 1*PS_devs_Fp_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_Fp_slow_SXL(:,3) - PS_range, 1*PS_devs_Fp_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3);
title('Slow (NXL fitting)','FontSize',8);
%xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-1.6 1.6]);

h2=subplot(4,4,5); % fast, FXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_fast_exclude(:,1) - PS_range, 1*PS_devs_Fp_fast_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_fast_exclude(:,2) - PS_range, 1*PS_devs_Fp_fast_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_Fp_fast_exclude(:,3) - PS_range, 1*PS_devs_Fp_fast_exclude(:,3),'LineWidth',1.1,'Color',Colour3);
ylabel({'{\bfw/ exclusion}'},'FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-1.6 1.6]);

subplot(4,4,6); % slow, FXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_slow_exclude(:,1) - PS_range, 1*PS_devs_Fp_slow_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_slow_exclude(:,2) - PS_range, 1*PS_devs_Fp_slow_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_Fp_slow_exclude(:,3) - PS_range, 1*PS_devs_Fp_slow_exclude(:,3),'LineWidth',1.1,'Color',Colour3);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-1.6 1.6]);

subplot(4,4,7); % fast, NXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_fast_SXL_exclude(:,1) - PS_range, 1*PS_devs_Fp_fast_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_fast_SXL_exclude(:,2) - PS_range, 1*PS_devs_Fp_fast_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_Fp_fast_SXL_exclude(:,3) - PS_range, 1*PS_devs_Fp_fast_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-1.6 1.6]);

subplot(4,4,8); % slow, NXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_slow_SXL_exclude(:,1) - PS_range, 1*PS_devs_Fp_slow_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_slow_SXL_exclude(:,2) - PS_range, 1*PS_devs_Fp_slow_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_Fp_slow_SXL_exclude(:,3) - PS_range, 1*PS_devs_Fp_slow_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-1.6 1.6]);

h3=subplot(4,4,9); % fast, FXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True v_p','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_Fp_fast(:,1) - vP_range, 1*vP_devs_Fp_fast(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_Fp_fast(:,2) - vP_range, 1*vP_devs_Fp_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_Fp_fast(:,3) - vP_range, 1*vP_devs_Fp_fast(:,3),'LineWidth',1.1,'Color',Colour3);
ylabel({'{\bfNo exclusion}'},'FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-1.3 0.2]);

subplot(4,4,10) % slow, FXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True v_p','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_Fp_slow(:,1) - vP_range, 1*vP_devs_Fp_slow(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_Fp_slow(:,2) - vP_range, 1*vP_devs_Fp_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_Fp_slow(:,3) - vP_range, 1*vP_devs_Fp_slow(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.6 0.2]);

subplot(4,4,11) % fast, NXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True v_p','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_Fp_fast_SXL(:,1) - vP_range, 1*vP_devs_Fp_fast_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_Fp_fast_SXL(:,2) - vP_range, 1*vP_devs_Fp_fast_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_Fp_fast_SXL(:,3) - vP_range, 1*vP_devs_Fp_fast_SXL(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.6 0.2]);

subplot(4,4,12) % slow, NXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True v_p','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_Fp_slow_SXL(:,1) - vP_range, 1*vP_devs_Fp_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_Fp_slow_SXL(:,2) - vP_range, 1*vP_devs_Fp_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_Fp_slow_SXL(:,3) - vP_range, 1*vP_devs_Fp_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.6 0.2]);

h4=subplot(4,4,13) % fast, FXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True v_p','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_Fp_fast_exclude(:,1) - vP_range, 1*vP_devs_Fp_fast_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_Fp_fast_exclude(:,2) - vP_range, 1*vP_devs_Fp_fast_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_Fp_fast_exclude(:,3) - vP_range, 1*vP_devs_Fp_fast_exclude(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.6 0.2]);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
ylabel({'{\bfw/ exclusion}'},'FontSize',8);

subplot(4,4,14) % slow, FXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True v_p','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_Fp_slow_exclude(:,1) - vP_range, 1*vP_devs_Fp_slow_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_Fp_slow_exclude(:,2) - vP_range, 1*vP_devs_Fp_slow_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_Fp_slow_exclude(:,3) - vP_range, 1*vP_devs_Fp_slow_exclude(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.6 0.2]);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);

subplot(4,4,15) % fast, NXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True v_p','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_Fp_fast_SXL_exclude(:,1) - vP_range, 1*vP_devs_Fp_fast_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_Fp_fast_SXL_exclude(:,2) - vP_range, 1*vP_devs_Fp_fast_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_Fp_fast_SXL_exclude(:,3) - vP_range, 1*vP_devs_Fp_fast_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.6 0.2]);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);

subplot(4,4,16) % slow, NXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True v_p','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_Fp_slow_SXL_exclude(:,1) - vP_range, 1*vP_devs_Fp_slow_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_Fp_slow_SXL_exclude(:,2) - vP_range, 1*vP_devs_Fp_slow_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_Fp_slow_SXL_exclude(:,3) - vP_range, 1*vP_devs_Fp_slow_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.6 0.2]);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);

p1=get(h1,'position');
p2=get(h2,'position');
height=p1(2)+p1(4)-p2(2);
hx1=axes('position',[0.11 p2(2) p2(3) height],'visible','off');
h_label=ylabel('Fitted {\itPS} error (x10^{-4} min^{-1} )','visible','on');

p3=get(h3,'position');
p4=get(h4,'position');
height=p3(2)+p3(4)-p4(2);
hx2=axes('position',[0.11 p4(2) p4(3) height],'visible','off');
h_label=ylabel('Fitted {\itv_p} error (x10^{-2})','visible','on');

set(gcf, 'units', 'centimeters','Position', [5 5 17.56 21.08]);

annotation(figure(1),'textbox',[0.090 0.918 0.05 0.045],'String','(A)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.297 0.918 0.06 0.045],'String','(B)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.503 0.918 0.06 0.045],'String','(C)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.709 0.918 0.06 0.045],'String','(D)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.090 0.697 0.06 0.045],'String','(E)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.297 0.697 0.06 0.045],'String','(F)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.503 0.697 0.06 0.045],'String','(G)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.709 0.697 0.06 0.045],'String','(H)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.090 0.476 0.06 0.045],'String','(I)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.297 0.476 0.06 0.045],'String','(J)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.503 0.476 0.06 0.045],'String','(K)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.709 0.476 0.06 0.045],'String','(L)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.090 0.255 0.06 0.045],'String','(M)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.297 0.255 0.06 0.045],'String','(N)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.503 0.255 0.06 0.045],'String','(O)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.709 0.255 0.06 0.045],'String','(P)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);

set(gcf, 'units', 'centimeters','PaperPosition', [0 0 17.56 21.08]);    % can be bigger than screen
print(gcf, 'Figure_3.png', '-dpng','-r800');
print(gcf, 'Figure_3.tif', '-dtiff','-r800');
