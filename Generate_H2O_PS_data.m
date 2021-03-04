
% Generates the water exchange PS simulation data as shown in Manning et al.
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;

addpath('DCE_Simulation_Functions');
    
% Select default parameters
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
 
% ranges of PS, vP, kbe to test
PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;
vP_range = linspace(SimParam.min_vP,SimParam.max_vP,10)';
PS_fixed = PhysParam.PS_fixed;

N_PS = size(PS_range,1); %range sizes to test
N_vP = size(vP_range,1); %range sizes to test
kbe_ranges = [0 1.375 2.75 5.5 1000]; % range of kbe values to test

%% Sim water exchange with Patlak fitting (fast injection, Patlak fit)
% T1 acquisition
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
    for i = 1:size(kbe_ranges,2);
        PhysParam.kbe_perS = kbe_ranges(i);
        for i_PS = 1:N_PS
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
            PhysParam.vP = vP_fixed(1);
            PhysParam.PS_perMin = PS_range(i_PS);
            [temp, PS_fit_2S1X_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        for i_vP = 1:N_vP
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
            PhysParam.PS = PS_fixed(1);
            PhysParam.vP = vP_range(i_vP);
            [vP_fit_2S1X_fast(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        
        PS_means_H2O_fast(:,i) = mean(PS_fit_2S1X_fast,1)'; % add mean for each PS 
        PS_devs_H2O_fast(:,i) = std(PS_fit_2S1X_fast,0,1)'; % add standard deviation
        
        vP_means_H2O_fast(:,i) = mean(vP_fit_2S1X_fast,1)'; % add mean for each vP
        vP_devs_H2O_fast(:,i) = std(vP_fit_2S1X_fast,0,1)'; % add standard deviation
    end
    
    %% Sim water exchange with Patlak fitting (fast injection, Patlak fit, w/ exclusion)
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
    for i = 1:size(kbe_ranges,2);
        PhysParam.kbe_perS = kbe_ranges(i);
        for i_PS = 1:N_PS
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
            PhysParam.vP = vP_fixed(1);
            PhysParam.PS_perMin = PS_range(i_PS);
            [temp, PS_fit_2S1X_fast_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        for i_vP = 1:N_vP
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
            PhysParam.PS = PS_fixed(1);
            PhysParam.vP = vP_range(i_vP);
            [vP_fit_2S1X_fast_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        
        PS_means_H2O_fast_exclude(:,i) = mean(PS_fit_2S1X_fast_exclude,1)'; % add mean for each PS 
        PS_devs_H2O_fast_exclude(:,i) = std(PS_fit_2S1X_fast_exclude,0,1)'; % add standard deviation
        
        vP_means_H2O_fast_exclude(:,i) = mean(vP_fit_2S1X_fast_exclude,1)'; % add mean for each vP
        vP_devs_H2O_fast_exclude(:,i) = std(vP_fit_2S1X_fast_exclude,0,1)'; % add standard deviation
    end
    %% Sim water exchange with Patlak fitting (fast injection, SXL fit)
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
    for i = 1:size(kbe_ranges,2);
        PhysParam.kbe_perS = kbe_ranges(i);
        for i_PS = 1:N_PS
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
            PhysParam.vP = vP_fixed(1);
            PhysParam.PS_perMin = PS_range(i_PS);
            [temp, PS_fit_2S1X_SXL(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        for i_vP = 1:N_vP
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
            PhysParam.PS = PS_fixed(1);
            PhysParam.vP = vP_range(i_vP);
            [vP_fit_2S1X_SXL(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        PS_means_H2O_SXL(:,i) = mean(PS_fit_2S1X_SXL,1)'; % add mean for each PS
        PS_devs_H2O_SXL(:,i) = std(PS_fit_2S1X_SXL,0,1)'; % add standard deviation
        
        vP_means_H2O_SXL(:,i) = mean(vP_fit_2S1X_SXL,1)'; % add mean for each vP
        vP_devs_H2O_SXL(:,i) = std(vP_fit_2S1X_SXL,0,1)'; % add standard deviation
    end
    
    %% Sim water exchange with Patlak fitting (fast injection, SXL fit, w/ exclusion)
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
    for i = 1:size(kbe_ranges,2);
        PhysParam.kbe_perS = kbe_ranges(i);
        for i_PS = 1:N_PS
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
            PhysParam.vP = vP_fixed(1);
            PhysParam.PS_perMin = PS_range(i_PS);
            [temp, PS_fit_2S1X_SXL_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        for i_vP = 1:N_vP
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
            PhysParam.PS = PS_fixed(1);
            PhysParam.vP = vP_range(i_vP);
            [vP_fit_2S1X_SXL_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        PS_means_H2O_SXL_exclude(:,i) = mean(PS_fit_2S1X_SXL_exclude,1)'; % add mean for each PS
        PS_devs_H2O_SXL_exclude(:,i) = std(PS_fit_2S1X_SXL_exclude,0,1)'; % add standard deviation
        
        vP_means_H2O_SXL_exclude(:,i) = mean(vP_fit_2S1X_SXL_exclude,1)'; % add mean for each vP
        vP_devs_H2O_SXL_exclude(:,i) = std(vP_fit_2S1X_SXL_exclude,0,1)'; % add standard deviation
    end

    %% Sim water exchange (slow injection, Patlak fit)
    SimParam.SXLfit = 0;
    SimParam.InjectionRate = 'slow';
    SimParam.t_start_s = 0;
    SimParam.NIgnore = max(SimParam.baselineScans);
    
    load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
    SimParam.Cp_AIF_mM = Cp_AIF_mM;
    SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
    SimParam.InputAIFDCENFrames = 32; % number of time points
    
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
    for i = 1:size(kbe_ranges,2);
        PhysParam.kbe_perS = kbe_ranges(i);
        for i_PS = 1:N_PS
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
            PhysParam.vP = vP_fixed(1);
            PhysParam.PS_perMin = PS_range(i_PS);
            [temp, PS_fit_2S1X_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        for i_vP = 1:N_vP
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
            PhysParam.PS = PS_fixed(1);
            PhysParam.vP = vP_range(i_vP);
            [vP_fit_2S1X_slow(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        PS_means_H2O_slow(:,i) = mean(PS_fit_2S1X_slow,1)'; % add mean for each PS
        PS_devs_H2O_slow(:,i) = std(PS_fit_2S1X_slow,0,1)'; % add standard deviation
        
        vP_means_H2O_slow(:,i) = mean(vP_fit_2S1X_slow,1)'; % add mean for each vP
        vP_devs_H2O_slow(:,i) = std(vP_fit_2S1X_slow,0,1)'; % add standard deviation
        
    end
    
    %% Sim water exchange (slow injection, Patlak fit, w/ exclusion)
    SimParam.SXLfit = 0;
    SimParam.InjectionRate = 'slow';
    SimParam.t_start_s = 0;
    SimParam.NIgnore = max(SimParam.baselineScans) + 3;
    
    load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
    SimParam.Cp_AIF_mM = Cp_AIF_mM;
    SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
    SimParam.InputAIFDCENFrames = 32; % number of time points

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
    for i = 1:size(kbe_ranges,2);
        PhysParam.kbe_perS = kbe_ranges(i);
        for i_PS = 1:N_PS
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
            PhysParam.vP = vP_fixed(1);
            PhysParam.PS_perMin = PS_range(i_PS);
            [temp, PS_fit_2S1X_slow_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        for i_vP = 1:N_vP
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
            PhysParam.PS = PS_fixed(1);
            PhysParam.vP = vP_range(i_vP);
            [vP_fit_2S1X_slow_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        PS_means_H2O_slow_exclude(:,i) = mean(PS_fit_2S1X_slow_exclude,1)'; % add mean for each PS
        PS_devs_H2O_slow_exclude(:,i) = std(PS_fit_2S1X_slow_exclude,0,1)'; % add standard deviation
        
        vP_means_H2O_slow_exclude(:,i) = mean(vP_fit_2S1X_slow_exclude,1)'; % add mean for each vP
        vP_devs_H2O_slow_exclude(:,i) = std(vP_fit_2S1X_slow_exclude,0,1)'; % add standard deviation
        
    end
    %% Sim water exchange (slow injection, SXL fit)
    SimParam.SXLfit = 1;
    
    for m = 1:N_PS
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
            [T1_blood_meas_s(n,1),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
            [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
    end
% loop through PS and vP values, simulate
    for i = 1:size(kbe_ranges,2);
        PhysParam.kbe_perS = kbe_ranges(i);
        for i_PS = 1:N_PS
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
            PhysParam.vP = vP_fixed(1);
            PhysParam.PS_perMin = PS_range(i_PS);
            [temp, PS_fit_2S1X_slow_SXL(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        for i_vP = 1:N_vP
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
            PhysParam.PS = PS_fixed(1);
            PhysParam.vP = vP_range(i_vP);
            [vP_fit_2S1X_slow(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        PS_means_H2O_slow_SXL(:,i) = mean(PS_fit_2S1X_slow_SXL,1)'; % add mean for each PS
        PS_devs_H2O_slow_SXL(:,i) = std(PS_fit_2S1X_slow_SXL,0,1)'; % add standard deviation
        
        vP_means_H2O_slow_SXL(:,i) = mean(vP_fit_2S1X_slow,1)'; % add mean for each vP
        vP_devs_H2O_slow_SXL(:,i) = std(vP_fit_2S1X_slow,0,1)'; % add standard deviation
    end

     %% Sim water exchange (slow injection, SXL fit, w/ exclusion)
    SimParam.SXLfit = 1;
    SimParam.NIgnore = max(SimParam.baselineScans) + 3;

    for m = 1:N_PS
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
            [T1_blood_meas_s(n,1),temp,acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
            [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
    end
% loop through PS and vP values, simulate
    for i = 1:size(kbe_ranges,2);
        PhysParam.kbe_perS = kbe_ranges(i);
        for i_PS = 1:N_PS
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
            PhysParam.vP = vP_fixed(1);
            PhysParam.PS_perMin = PS_range(i_PS);
            [temp, PS_fit_2S1X_slow_SXL_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        for i_vP = 1:N_vP
            PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
            PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
            PhysParam.PS = PS_fixed(1);
            PhysParam.vP = vP_range(i_vP);
            [vP_fit_2S1X_slow_exclude(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
        end
        PS_means_H2O_slow_SXL_exclude(:,i) = mean(PS_fit_2S1X_slow_SXL_exclude,1)'; % add mean for each PS
        PS_devs_H2O_slow_SXL_exclude(:,i) = std(PS_fit_2S1X_slow_SXL_exclude,0,1)'; % add standard deviation
        
        vP_means_H2O_slow_SXL_exclude(:,i) = mean(vP_fit_2S1X_slow_exclude,1)'; % add mean for each vP
        vP_devs_H2O_slow_SXL_exclude(:,i) = std(vP_fit_2S1X_slow_exclude,0,1)'; % add standard deviation
    end

    %% Save simulation data
    PS_range = PS_range * 1e4;
    PS_means_H2O_fast = PS_means_H2O_fast * 1e4;
    PS_means_H2O_fast_exclude = PS_means_H2O_fast_exclude * 1e4;
    PS_means_H2O_SXL = PS_means_H2O_SXL * 1e4;
    PS_means_H2O_SXL_exclude = PS_means_H2O_SXL_exclude * 1e4;
    PS_means_H2O_slow = PS_means_H2O_slow * 1e4;
    PS_means_H2O_slow_exclude = PS_means_H2O_slow_exclude * 1e4;
    PS_means_H2O_slow_SXL = PS_means_H2O_slow_SXL * 1e4;
    PS_means_H2O_slow_SXL_exclude = PS_means_H2O_slow_SXL_exclude * 1e4;
    PS_devs_H2O_fast = PS_devs_H2O_fast * 1e4;
    PS_devs_H2O_fast_exclude = PS_devs_H2O_fast_exclude * 1e4;
    PS_devs_H2O_SXL = PS_devs_H2O_SXL * 1e4;
    PS_devs_H2O_SXL_exclude = PS_devs_H2O_SXL_exclude * 1e4;
    PS_devs_H2O_slow = PS_devs_H2O_slow * 1e4;
    PS_devs_H2O_slow_exclude = PS_devs_H2O_slow_exclude * 1e4;
    PS_devs_H2O_slow_SXL = PS_devs_H2O_slow_SXL * 1e4;
    PS_devs_H2O_slow_SXL_exclude = PS_devs_H2O_slow_SXL_exclude * 1e4;
    save('PS_means_H2O','PS_means_H2O_fast','PS_means_H2O_fast_exclude','PS_means_H2O_SXL',...
        'PS_means_H2O_SXL_exclude','PS_means_H2O_slow','PS_means_H2O_slow_exclude','PS_means_H2O_slow_SXL','PS_means_H2O_slow_SXL_exclude');
    save('PS_devs_H2O','PS_devs_H2O_fast','PS_devs_H2O_fast_exclude','PS_devs_H2O_SXL',...
        'PS_devs_H2O_SXL_exclude','PS_devs_H2O_slow','PS_devs_H2O_slow_exclude','PS_devs_H2O_slow_SXL','PS_devs_H2O_slow_SXL_exclude');
    
    vP_range = vP_range * 1e2;
    vP_means_H2O_fast = vP_means_H2O_fast * 1e2;
    vP_means_H2O_fast_exclude = vP_means_H2O_fast_exclude * 1e2;
    vP_means_H2O_SXL = vP_means_H2O_SXL * 1e2;
    vP_means_H2O_SXL_exclude = vP_means_H2O_SXL_exclude * 1e2;
    vP_means_H2O_slow = vP_means_H2O_slow * 1e2;
    vP_means_H2O_slow_exclude = vP_means_H2O_slow_exclude * 1e2;
    vP_means_H2O_slow_SXL = vP_means_H2O_slow_SXL * 1e2;
    vP_means_H2O_slow_SXL_exclude = vP_means_H2O_slow_SXL_exclude * 1e2;
    vP_devs_H2O_fast = vP_devs_H2O_fast * 1e2;
    vP_devs_H2O_fast_exclude = vP_devs_H2O_fast_exclude * 1e2;
    vP_devs_H2O_SXL = vP_devs_H2O_SXL * 1e2;
    vP_devs_H2O_SXL_exclude = vP_devs_H2O_SXL_exclude * 1e2;
    vP_devs_H2O_slow = vP_devs_H2O_slow * 1e2;
    vP_devs_H2O_slow_exclude = vP_devs_H2O_slow_exclude * 1e2;
    vP_devs_H2O_slow_SXL = vP_devs_H2O_slow_SXL * 1e2;
    vP_devs_H2O_slow_SXL_exclude = vP_devs_H2O_slow_SXL_exclude * 1e2;
    save('vP_means_H2O','vP_means_H2O_fast','vP_means_H2O_fast_exclude','vP_means_H2O_SXL',...
        'vP_means_H2O_SXL_exclude','vP_means_H2O_slow','vP_means_H2O_slow_exclude','vP_means_H2O_slow_SXL','vP_means_H2O_slow_SXL_exclude');
    save('vP_devs_H2O','vP_devs_H2O_fast','vP_devs_H2O_fast_exclude','vP_devs_H2O_SXL',...
       'vP_devs_H2O_SXL_exclude','vP_devs_H2O_slow','vP_devs_H2O_slow_exclude','vP_devs_H2O_slow_SXL','vP_devs_H2O_slow_SXL_exclude');

%% Plot figures of PS     
Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];
Colour4 = [0.494 0.184 0.556 0.5];
Colour5 = [0.466 0.674 0.188 0.5];

figure(1)

h1=subplot(2,4,1) % fast, FXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_fast(:,1) - PS_range, 1*PS_devs_H2O_fast(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(PS_range + 0.04, PS_means_H2O_fast(:,2) - PS_range, 1*PS_devs_H2O_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.08, PS_means_H2O_fast(:,3) - PS_range, 1*PS_devs_H2O_fast(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_fast(:,4) - PS_range, 1*PS_devs_H2O_fast(:,4),'Linewidth',1.1,'Color',Colour4); hold on;
errorbar(PS_range + 0.16, PS_means_H2O_fast(:,5) - PS_range, 1*PS_devs_H2O_fast(:,5),'Linewidth',1.1,'Color',Colour5, 'Linestyle', '--');
ylabel({'{\bfNo exclusion}'},'FontSize',8);
title('Bolus (FXL fitting)','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-5 5]);
legend({'{\itk_{be}} = 0 s^{-1}','{\itk_{be}} = 1.375 s^{-1}','{\itk_{be}} = 2.75 s^{-1}','{\itk_{be}} = 5.5 s^{-1}','{\itk_{be}} = 1000 s^{-1}'},'Position',[0.172 0.660 0.081 0.0207],'FontSize',8)
legend('boxoff')

subplot(2,4,2) % slow, FXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_slow(:,1) - PS_range, 1*PS_devs_H2O_slow(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(PS_range + 0.04, PS_means_H2O_slow(:,2) - PS_range, 1*PS_devs_H2O_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.08, PS_means_H2O_slow(:,3) - PS_range, 1*PS_devs_H2O_slow(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_slow(:,4) - PS_range, 1*PS_devs_H2O_slow(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(PS_range + 0.16, PS_means_H2O_slow(:,5) - PS_range, 1*PS_devs_H2O_slow(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
title('Slow inj (FXL fitting)','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2.75 2.75]);

subplot(2,4,3) % fast, NXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_SXL(:,1) - PS_range, 1*PS_devs_H2O_SXL(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(PS_range + 0.04, PS_means_H2O_SXL(:,2) - PS_range, 1*PS_devs_H2O_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.08, PS_means_H2O_SXL(:,3) - PS_range, 1*PS_devs_H2O_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_SXL(:,4) - PS_range, 1*PS_devs_H2O_SXL(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(PS_range + 0.16, PS_means_H2O_SXL(:,5) - PS_range, 1*PS_devs_H2O_SXL(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
title('Bolus (NXL fitting)','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2.75 2.75]);

subplot(2,4,4) % slow, NXL

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_slow_SXL(:,1) - PS_range, 1*PS_devs_H2O_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(PS_range + 0.04, PS_means_H2O_slow_SXL(:,2) - PS_range, 1*PS_devs_H2O_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.08, PS_means_H2O_slow_SXL(:,3) - PS_range, 1*PS_devs_H2O_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_slow_SXL(:,4) - PS_range, 1*PS_devs_H2O_slow_SXL(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(PS_range + 0.16, PS_means_H2O_slow_SXL(:,5) - PS_range, 1*PS_devs_H2O_slow_SXL(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
title('Slow inj (NXL fitting)','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2.75 2.75]);

h2=subplot(2,4,5) %, fast, FXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_fast_exclude(:,1) - PS_range, 1*PS_devs_H2O_fast_exclude(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(PS_range + 0.04, PS_means_H2O_fast_exclude(:,2) - PS_range, 1*PS_devs_H2O_fast_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.08, PS_means_H2O_fast_exclude(:,3) - PS_range, 1*PS_devs_H2O_fast_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_fast_exclude(:,4) - PS_range, 1*PS_devs_H2O_fast_exclude(:,4),'Linewidth',1.1,'Color',Colour4); hold on;
errorbar(PS_range + 0.16, PS_means_H2O_fast_exclude(:,5) - PS_range, 1*PS_devs_H2O_fast_exclude(:,5),'Linewidth',1.1,'Color',Colour5, 'Linestyle', '--');
ylabel({'{\bfw/ exclusion}'},'FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2.75 2.75]);

subplot(2,4,6) % slow, FXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_slow_exclude(:,1) - PS_range, 1*PS_devs_H2O_slow_exclude(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(PS_range + 0.04, PS_means_H2O_slow_exclude(:,2) - PS_range, 1*PS_devs_H2O_slow_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.08, PS_means_H2O_slow_exclude(:,3) - PS_range, 1*PS_devs_H2O_slow_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_slow_exclude(:,4) - PS_range, 1*PS_devs_H2O_slow_exclude(:,4),'Linewidth',1.1,'Color',Colour4); hold on;
errorbar(PS_range + 0.16, PS_means_H2O_slow_exclude(:,5) - PS_range, 1*PS_devs_H2O_slow_exclude(:,5),'Linewidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2.75 2.75]);

subplot(2,4,7) % fast, NXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_SXL_exclude(:,1) - PS_range, 1*PS_devs_H2O_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(PS_range + 0.04, PS_means_H2O_SXL_exclude(:,2) - PS_range, 1*PS_devs_H2O_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.08, PS_means_H2O_SXL_exclude(:,3) - PS_range, 1*PS_devs_H2O_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_SXL_exclude(:,4) - PS_range, 1*PS_devs_H2O_SXL_exclude(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(PS_range + 0.16, PS_means_H2O_SXL_exclude(:,5) - PS_range, 1*PS_devs_H2O_SXL_exclude(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlim([0 max(PS_range)+0.12]);
ylim([-2.75 2.75]);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);

subplot(2,4,8) % slow, NXL, exclude

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_slow_SXL_exclude(:,1) - PS_range, 1*PS_devs_H2O_slow_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(PS_range + 0.04, PS_means_H2O_slow_SXL_exclude(:,2) - PS_range, 1*PS_devs_H2O_slow_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.08, PS_means_H2O_slow_SXL_exclude(:,3) - PS_range, 1*PS_devs_H2O_slow_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_slow_SXL_exclude(:,4) - PS_range, 1*PS_devs_H2O_slow_SXL_exclude(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(PS_range + 0.16, PS_means_H2O_slow_SXL_exclude(:,5) - PS_range, 1*PS_devs_H2O_slow_SXL_exclude(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2.75 2.75]);

figure(1)
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
set(gcf,'PaperPositionMode','manual');
print(gcf, 'Figure_2.png', '-dpng','-r800');
print(gcf, 'Figure_2.tif', '-dtiff','-r800');

%% Plot Figure of vp (for supplementary material)
figure(2)
h3=subplot(2,4,1) % fast, FXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_fast(:,1) - vP_range, 1*vP_devs_H2O_fast(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(vP_range + 0.011, vP_means_H2O_fast(:,2) - vP_range, 1*vP_devs_H2O_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.022, vP_means_H2O_fast(:,3) - vP_range, 1*vP_devs_H2O_fast(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.033, vP_means_H2O_fast(:,4) - vP_range, 1*vP_devs_H2O_fast(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(vP_range + 0.044, vP_means_H2O_fast(:,5) - vP_range, 1*vP_devs_H2O_fast(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
ylabel({'{\bfNo exclusion}'},'FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-1.4 0.2]);
title('Bolus (FXL fitting)','FontSize',8);

subplot(2,4,2) % slow, FXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_slow(:,1) - vP_range, 1*vP_devs_H2O_slow(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(vP_range + 0.011, vP_means_H2O_slow(:,2) - vP_range, 1*vP_devs_H2O_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.022, vP_means_H2O_slow(:,3) - vP_range, 1*vP_devs_H2O_slow(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.033, vP_means_H2O_slow(:,4) - vP_range, 1*vP_devs_H2O_slow(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(vP_range + 0.044, vP_means_H2O_slow(:,5) - vP_range, 1*vP_devs_H2O_slow(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.8 0.8]);
title('Slow inj (FXL fitting)','FontSize',8);

subplot(2,4,3) % fast, NXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_SXL(:,1) - vP_range, 1*vP_devs_H2O_SXL(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(vP_range + 0.011, vP_means_H2O_SXL(:,2) - vP_range, 1*vP_devs_H2O_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.022, vP_means_H2O_SXL(:,3) - vP_range, 1*vP_devs_H2O_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.033, vP_means_H2O_SXL(:,4) - vP_range, 1*vP_devs_H2O_SXL(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(vP_range + 0.044, vP_means_H2O_SXL(:,5) - vP_range, 1*vP_devs_H2O_SXL(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.8 0.8]);
title('Bolus (NXL fitting)','FontSize',8);

subplot(2,4,4) % slow, NXL

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_slow_SXL(:,1) - vP_range, 1*vP_devs_H2O_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(vP_range + 0.011, vP_means_H2O_slow_SXL(:,2) - vP_range, 1*vP_devs_H2O_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.022, vP_means_H2O_slow_SXL(:,3) - vP_range, 1*vP_devs_H2O_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.033, vP_means_H2O_slow_SXL(:,4) - vP_range, 1*vP_devs_H2O_slow_SXL(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(vP_range + 0.044, vP_means_H2O_slow_SXL(:,5) - vP_range, 1*vP_devs_H2O_slow_SXL(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.8 0.8]);
title('Slow inj (NXL fitting)','FontSize',8);

h4=subplot(2,4,5) % fast, FXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_fast_exclude(:,1) - vP_range, 1*vP_devs_H2O_fast_exclude(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(vP_range + 0.011, vP_means_H2O_fast_exclude(:,2) - vP_range, 1*vP_devs_H2O_fast_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.022, vP_means_H2O_fast_exclude(:,3) - vP_range, 1*vP_devs_H2O_fast_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.033, vP_means_H2O_fast_exclude(:,4) - vP_range, 1*vP_devs_H2O_fast_exclude(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(vP_range + 0.044, vP_means_H2O_fast_exclude(:,5) - vP_range, 1*vP_devs_H2O_fast_exclude(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
ylabel({'{\bfw/ exclusion}'},'FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.8 0.8]);
legend({'{\itk_{be}} = 0 s^{-1}','{\itk_{be}} = 1.375 s^{-1}','{\itk_{be}} = 2.75 s^{-1}','{\itk_{be}} = 5.5 s^{-1}','{\itk_{be}} = 1000 s^{-1}'},'Position',[0.172 0.355 0.081 0.0207],'FontSize',8)
legend('boxoff')

subplot(2,4,6) % slow, FXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_slow_exclude(:,1) - vP_range, 1*vP_devs_H2O_slow_exclude(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(vP_range + 0.011, vP_means_H2O_slow_exclude(:,2) - vP_range, 1*vP_devs_H2O_slow_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.022, vP_means_H2O_slow_exclude(:,3) - vP_range, 1*vP_devs_H2O_slow_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.033, vP_means_H2O_slow_exclude(:,4) - vP_range, 1*vP_devs_H2O_slow_exclude(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(vP_range + 0.044, vP_means_H2O_slow_exclude(:,5) - vP_range, 1*vP_devs_H2O_slow_exclude(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.8 0.8]);

subplot(2,4,7) % fast, NXL, exclude

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_SXL_exclude(:,1) - vP_range, 1*vP_devs_H2O_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(vP_range + 0.011, vP_means_H2O_SXL_exclude(:,2) - vP_range, 1*vP_devs_H2O_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.022, vP_means_H2O_SXL_exclude(:,3) - vP_range, 1*vP_devs_H2O_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.033, vP_means_H2O_SXL_exclude(:,4) - vP_range, 1*vP_devs_H2O_SXL_exclude(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(vP_range + 0.044, vP_means_H2O_SXL_exclude(:,5) - vP_range, 1*vP_devs_H2O_SXL_exclude(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlim([min(vP_range) max(vP_range)+0.026]);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
ylim([-0.8 0.8]);

subplot(2,4,8) % slow, NXL, exclude
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_slow_SXL_exclude(:,1) - vP_range, 1*vP_devs_H2O_slow_SXL_exclude(:,1),'LineWidth',1.1,'Color',Colour1, 'Linestyle', '--'); hold on;
errorbar(vP_range + 0.011, vP_means_H2O_slow_SXL_exclude(:,2) - vP_range, 1*vP_devs_H2O_slow_SXL_exclude(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.022, vP_means_H2O_slow_SXL_exclude(:,3) - vP_range, 1*vP_devs_H2O_slow_SXL_exclude(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
errorbar(vP_range + 0.033, vP_means_H2O_slow_SXL_exclude(:,4) - vP_range, 1*vP_devs_H2O_slow_SXL_exclude(:,4),'LineWidth',1.1,'Color',Colour4); hold on;
errorbar(vP_range + 0.044, vP_means_H2O_slow_SXL_exclude(:,5) - vP_range, 1*vP_devs_H2O_slow_SXL_exclude(:,5),'LineWidth',1.1,'Color',Colour5, 'Linestyle', '--');
xlim([min(vP_range) max(vP_range)+0.026]);
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
ylim([-0.8 0.8]);
set(gcf, 'units', 'centimeters','Position', [5 5 20 16]);

figure(2)
p3=get(h3,'position');
p4=get(h4,'position');
height=p3(2)+p3(4)-p4(2);
hx2=axes('position',[0.11 p4(2) p4(3) height],'visible','off');
h_label=ylabel('Fitted {\itv_p} error (x10^{-2})','visible','on');

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
print(gcf, 'Supp_Figure_vP_H2O.png', '-dpng','-r800');
