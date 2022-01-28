% Generates the B1 inhomogeneity PS simulation data as shown in Manning et al.
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;

addpath('DCE_Simulation_Functions');

[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;

PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;
vP_range = linspace(SimParam.min_vP,SimParam.max_vP,10)';
PS_fixed = PhysParam.PS_fixed;

blood_k_values = [1 1.12 0.95]; % blood K_FA values
tissue_k_values = [1 1.12 1.12]; % tissue K_FA values
N_PS = size(PS_range,1); %range sizes to test
N_vP = size(vP_range,1); %range sizes to test

%% Generate B1 inhomogeneity sims
%% B1 inhomogeneity figures (fast injection, Patlak fitting)
SimParam.SXLfit = 0; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;

for k = 1:size(blood_k_values,2)
    DCESeqParam.blood_FA_error = blood_k_values(k);
    DCESeqParam.tissue_FA_error = tissue_k_values(k);
    T1acqParam.blood_FA_true_rads = DCESeqParam.blood_FA_error * T1acqParam.blood_FA_nom_rads; % Actual FA experienced by blood in T1 acquisition
    T1acqParam.tissue_FA_true_rads = DCESeqParam.tissue_FA_error * T1acqParam.tissue_FA_nom_rads; % Actual FA experienced by tissue in T1 acquisition
    DCESeqParam.blood_FA_true_deg = DCESeqParam.blood_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by blood in DCE
    DCESeqParam.tissue_FA_true_deg = DCESeqParam.tissue_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by tissue in DCE
    DCESeqParam.blood_FA_meas_deg = DCESeqParam.FA_nom_deg;
    DCESeqParam.tissue_FA_meas_deg = DCESeqParam.FA_nom_deg;
    T1acqParam.T1_SNR = 318;
    
    % loop through PS and vP values, simulate
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
            [T1_blood_meas_s(n,i_PS),temp,acqParam.FA_error_meas(n,i_PS),temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    for i_vP = 1:N_vP
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
        PhysParam.PS = PS_fixed(1);
        PhysParam.vP = vP_range(i_vP);
        [vP_fit_fast(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_B1_fast(:,k) = mean(PS_fit_fast,1)'; % mean for each PS
    PS_devs_B1_fast(:,k) = std(PS_fit_fast,0,1)'; % standard deviation for each PS
    
    vP_means_B1_fast(:,k) = mean(vP_fit_fast,1)'; % mean for each vP
    vP_devs_B1_fast(:,k) = std(vP_fit_fast,0,1)'; % standard deviation for each vP
end
%% B1 inhomogeneity figures (fast injection, Patlak fitting, B1 corrected)
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.SXLfit = 0; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;

for k = 1:size(blood_k_values,2)
    DCESeqParam.blood_FA_error = blood_k_values(k);
    DCESeqParam.tissue_FA_error = tissue_k_values(k);
    T1acqParam.blood_FA_true_rads = DCESeqParam.blood_FA_error * T1acqParam.blood_FA_nom_rads; % Actual FA experienced by blood in T1 acquisition
    T1acqParam.tissue_FA_true_rads = DCESeqParam.tissue_FA_error * T1acqParam.tissue_FA_nom_rads; % Actual FA experienced by tissue in T1 acquisition
    DCESeqParam.blood_FA_true_deg = DCESeqParam.blood_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by blood in DCE
    DCESeqParam.tissue_FA_true_deg = DCESeqParam.tissue_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by tissue in DCE
    DCESeqParam.blood_FA_meas_deg = DCESeqParam.blood_FA_true_deg; % If B1 correction is on (for DCE) correct FA
    DCESeqParam.tissue_FA_meas_deg = DCESeqParam.tissue_FA_true_deg; % If B1 correction is on (for DCE) correct FA
    T1acqParam.T1_SNR = 318;
    
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_true_rads;
            [T1_blood_meas_s(n,i_PS),temp,acqParam.FA_error_meas(n,i_PS),temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_true_rads;
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_fast_corrected(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_fast_corrected(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_B1_fast_corrected(:,k) = mean(PS_fit_fast_corrected,1)'; % mean for each PS
    PS_devs_B1_fast_corrected(:,k) = std(PS_fit_fast_corrected,0,1)'; % standard deviation for each PS
    
    vP_means_B1_fast_corrected(:,k) = mean(vP_fit_fast_corrected,1)'; % mean for each vP
    vP_devs_B1_fast_corrected(:,k) = std(vP_fit_fast_corrected,0,1)'; % standard deviation for each vP
end

%% B1 inhomogeneity figures (fast injection, NXL fitting)
SimParam.SXLfit = 1; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;

for k = 1:size(blood_k_values,2)
    DCESeqParam.blood_FA_error = blood_k_values(k);
    DCESeqParam.tissue_FA_error = tissue_k_values(k);
    T1acqParam.blood_FA_true_rads = DCESeqParam.blood_FA_error * T1acqParam.blood_FA_nom_rads; % Actual FA experienced by blood in T1 acquisition
    T1acqParam.tissue_FA_true_rads = DCESeqParam.tissue_FA_error * T1acqParam.tissue_FA_nom_rads; % Actual FA experienced by tissue in T1 acquisition
    DCESeqParam.blood_FA_true_deg = DCESeqParam.blood_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by blood in DCE
    DCESeqParam.tissue_FA_true_deg = DCESeqParam.tissue_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by tissue in DCE
    DCESeqParam.blood_FA_meas_deg = DCESeqParam.FA_nom_deg;
    DCESeqParam.tissue_FA_meas_deg = DCESeqParam.FA_nom_deg;
    T1acqParam.T1_SNR = 318;
    
    % loop through PS and vP values, simulate
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
            [T1_blood_meas_s(n,i_PS),temp,acqParam.FA_error_meas(n,i_PS),temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_fast_SXL(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    for i_vP = 1:N_vP
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
        PhysParam.PS = PS_fixed(1);
        PhysParam.vP = vP_range(i_vP);
        [vP_fit_fast_SXL(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_B1_fast_SXL(:,k) = mean(PS_fit_fast_SXL,1)'; % mean for each PS
    PS_devs_B1_fast_SXL(:,k) = std(PS_fit_fast_SXL,0,1)'; % standard deviation for each PS
    
    vP_means_B1_fast_SXL(:,k) = mean(vP_fit_fast_SXL,1)'; % mean for each vP
    vP_devs_B1_fast_SXL(:,k) = std(vP_fit_fast_SXL,0,1)'; % standard deviation for each vP
end

%% B1 inhomogeneity figures (fast injection, NXL fitting, B1 corrected)
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.SXLfit = 1; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;

for k = 1:size(blood_k_values,2)
    DCESeqParam.blood_FA_error = blood_k_values(k);
    DCESeqParam.tissue_FA_error = tissue_k_values(k);
    T1acqParam.blood_FA_true_rads = DCESeqParam.blood_FA_error * T1acqParam.blood_FA_nom_rads; % Actual FA experienced by blood in T1 acquisition
    T1acqParam.tissue_FA_true_rads = DCESeqParam.tissue_FA_error * T1acqParam.tissue_FA_nom_rads; % Actual FA experienced by tissue in T1 acquisition
    DCESeqParam.blood_FA_true_deg = DCESeqParam.blood_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by blood in DCE
    DCESeqParam.tissue_FA_true_deg = DCESeqParam.tissue_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by tissue in DCE
    DCESeqParam.blood_FA_meas_deg = DCESeqParam.blood_FA_true_deg; % If B1 correction is on (for DCE) correct FA
    DCESeqParam.tissue_FA_meas_deg = DCESeqParam.tissue_FA_true_deg; % If B1 correction is on (for DCE) correct FA
    T1acqParam.T1_SNR = 318;
    
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_true_rads;
            [T1_blood_meas_s(n,i_PS),temp,acqParam.FA_error_meas(n,i_PS),temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_true_rads;
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_fast_SXL_corrected(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_fast_SXL_corrected(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_B1_fast_SXL_corrected(:,k) = mean(PS_fit_fast_SXL_corrected,1)'; % mean for each PS
    PS_devs_B1_fast_SXL_corrected(:,k) = std(PS_fit_fast_SXL_corrected,0,1)'; % standard deviation for each PS
    
    vP_means_B1_fast_SXL_corrected(:,k) = mean(vP_fit_fast_SXL_corrected,1)'; % mean for each vP
    vP_devs_B1_fast_SXL_corrected(:,k) = std(vP_fit_fast_SXL_corrected,0,1)'; % standard deviation for each vP
end

%% B1 inhomogeneity figures (slow injection)
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.SXLfit = 0; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 SimParam.t_start_s = 0;
 SimParam.InjectionRate = 'slow';
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
 
 for k = 1:size(blood_k_values,2)
    DCESeqParam.blood_FA_error = blood_k_values(k);
    DCESeqParam.tissue_FA_error = tissue_k_values(k);
    T1acqParam.blood_FA_true_rads = DCESeqParam.blood_FA_error * T1acqParam.blood_FA_nom_rads; % Actual FA experienced by blood in T1 acquisition
    T1acqParam.tissue_FA_true_rads = DCESeqParam.tissue_FA_error * T1acqParam.tissue_FA_nom_rads; % Actual FA experienced by tissue in T1 acquisition
    DCESeqParam.blood_FA_true_deg = DCESeqParam.blood_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by blood in DCE
    DCESeqParam.tissue_FA_true_deg = DCESeqParam.tissue_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by tissue in DCE
    DCESeqParam.blood_FA_meas_deg = DCESeqParam.FA_nom_deg;
    DCESeqParam.tissue_FA_meas_deg = DCESeqParam.FA_nom_deg;
    T1acqParam.T1_SNR = 318;
    
    % loop through PS and vP values, simulate
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
            [T1_blood_meas_s(n,i_PS),temp,acqParam.FA_error_meas(n,i_PS),temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_slow(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
     
     PS_means_B1_slow(:,k) = mean(PS_fit_slow,1)'; % mean for each PS
     PS_devs_B1_slow(:,k) = std(PS_fit_slow,0,1)'; % standard deviation of PS
     
    vP_means_B1_slow(:,k) = mean(vP_fit_slow,1)'; % mean for each vP
    vP_devs_B1_slow(:,k) = std(vP_fit_slow,0,1)'; % standard deviation for each vP
 end

%% Generate B1 inhomogeneity figures (slow injection - B1 corrected)
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.SXLfit = 0; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 SimParam.t_start_s = 0;
 SimParam.InjectionRate = 'slow';
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
for k = 1:size(blood_k_values,2)
    DCESeqParam.blood_FA_error = blood_k_values(k);
    DCESeqParam.tissue_FA_error = tissue_k_values(k);
    T1acqParam.blood_FA_true_rads = DCESeqParam.blood_FA_error * T1acqParam.blood_FA_nom_rads; % Actual FA experienced by blood in T1 acquisition
    T1acqParam.tissue_FA_true_rads = DCESeqParam.tissue_FA_error * T1acqParam.tissue_FA_nom_rads; % Actual FA experienced by tissue in T1 acquisition
    DCESeqParam.blood_FA_true_deg = DCESeqParam.blood_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by blood in DCE
    DCESeqParam.tissue_FA_true_deg = DCESeqParam.tissue_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by tissue in DCE
    DCESeqParam.blood_FA_meas_deg = DCESeqParam.blood_FA_true_deg; % If B1 correction is on (for DCE) correct FA
    DCESeqParam.tissue_FA_meas_deg = DCESeqParam.tissue_FA_true_deg; % If B1 correction is on (for DCE) correct FA
    T1acqParam.T1_SNR = 318;
    
    % loop through PS and vP values, simulate
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_true_rads;
            [T1_blood_meas_s(n,i_PS),temp,acqParam.FA_error_meas(n,i_PS),temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_true_rads;
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_slow_corrected(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);    
    end
    for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_slow_corrected(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_B1_slow_corrected(:,k) = mean(PS_fit_slow_corrected,1)'; % mean for each PS
    PS_devs_B1_slow_corrected(:,k) = std(PS_fit_slow_corrected,0,1)'; % standard deviation for each PS
    
    vP_means_B1_slow_corrected(:,k) = mean(vP_fit_slow_corrected,1)'; % mean for each vP
    vP_devs_B1_slow_corrected(:,k) = std(vP_fit_slow_corrected,0,1)'; % standard deviation for each vP
end

%% B1 inhomogeneity figures (slow injection, NXL, no correction)
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.SXLfit = 1; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 SimParam.t_start_s = 0;
 SimParam.InjectionRate = 'slow';
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
 
 for k = 1:size(blood_k_values,2)
    DCESeqParam.blood_FA_error = blood_k_values(k);
    DCESeqParam.tissue_FA_error = tissue_k_values(k);
    T1acqParam.blood_FA_true_rads = DCESeqParam.blood_FA_error * T1acqParam.blood_FA_nom_rads; % Actual FA experienced by blood in T1 acquisition
    T1acqParam.tissue_FA_true_rads = DCESeqParam.tissue_FA_error * T1acqParam.tissue_FA_nom_rads; % Actual FA experienced by tissue in T1 acquisition
    DCESeqParam.blood_FA_true_deg = DCESeqParam.blood_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by blood in DCE
    DCESeqParam.tissue_FA_true_deg = DCESeqParam.tissue_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by tissue in DCE
    DCESeqParam.blood_FA_meas_deg = DCESeqParam.FA_nom_deg;
    DCESeqParam.tissue_FA_meas_deg = DCESeqParam.FA_nom_deg;
    T1acqParam.T1_SNR = 318;
    
    % loop through PS and vP values, simulate
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_nom_rads;
            [T1_blood_meas_s(n,i_PS),temp,acqParam.FA_error_meas(n,i_PS),temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_nom_rads;
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_slow_SXL(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_vP);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_slow_SXL(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
     
     PS_means_B1_slow_SXL(:,k) = mean(PS_fit_slow_SXL,1)'; % mean for each PS
     PS_devs_B1_slow_SXL(:,k) = std(PS_fit_slow_SXL,0,1)'; % standard deviation of PS
     
    vP_means_B1_slow_SXL(:,k) = mean(vP_fit_slow_SXL,1)'; % mean for each vP
    vP_devs_B1_slow_SXL(:,k) = std(vP_fit_slow_SXL,0,1)'; % standard deviation for each vP
 end
 
 %% Generate B1 inhomogeneity figures (slow injection - B1 corrected)
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.SXLfit = 1; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 SimParam.t_start_s = 0;
 SimParam.InjectionRate = 'slow';
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
for k = 1:size(blood_k_values,2)
    DCESeqParam.blood_FA_error = blood_k_values(k);
    DCESeqParam.tissue_FA_error = tissue_k_values(k);
    T1acqParam.blood_FA_true_rads = DCESeqParam.blood_FA_error * T1acqParam.blood_FA_nom_rads; % Actual FA experienced by blood in T1 acquisition
    T1acqParam.tissue_FA_true_rads = DCESeqParam.tissue_FA_error * T1acqParam.tissue_FA_nom_rads; % Actual FA experienced by tissue in T1 acquisition
    DCESeqParam.blood_FA_true_deg = DCESeqParam.blood_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by blood in DCE
    DCESeqParam.tissue_FA_true_deg = DCESeqParam.tissue_FA_error * DCESeqParam.FA_nom_deg; % Actual FA experienced by tissue in DCE
    DCESeqParam.blood_FA_meas_deg = DCESeqParam.blood_FA_true_deg; % If B1 correction is on (for DCE) correct FA
    DCESeqParam.tissue_FA_meas_deg = DCESeqParam.tissue_FA_true_deg; % If B1 correction is on (for DCE) correct FA
    T1acqParam.T1_SNR = 318;
    
    % loop through PS and vP values, simulate
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            T1acqParam.FA_true_rads = T1acqParam.blood_FA_true_rads;  % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.blood_FA_true_rads;
            [T1_blood_meas_s(n,i_PS),temp,acqParam.FA_error_meas(n,i_PS),temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
            T1acqParam.FA_true_rads = T1acqParam.tissue_FA_true_rads; % Seperate FA_true and FA_nom for blood and tissue
            T1acqParam.FA_nom_rads = T1acqParam.tissue_FA_true_rads;
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
        PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_slow_SXL_corrected(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);    
    end
    for i_vP = 1:N_vP
         PhysParam.T1_blood_meas_s = T1_blood_meas_s(:,i_PS);
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_slow_SXL_corrected(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_B1_slow_SXL_corrected(:,k) = mean(PS_fit_slow_SXL_corrected,1)'; % mean for each PS
    PS_devs_B1_slow_SXL_corrected(:,k) = std(PS_fit_slow_SXL_corrected,0,1)'; % standard deviation for each PS
    
    vP_means_B1_slow_SXL_corrected(:,k) = mean(vP_fit_slow_SXL_corrected,1)'; % mean for each vP
    vP_devs_B1_slow_SXL_corrected(:,k) = std(vP_fit_slow_SXL_corrected,0,1)'; % standard deviation for each vP
end

%% Graphs and scales
PS_range = PS_range * 1e4;
PS_means_B1_fast = PS_means_B1_fast * 1e4;
PS_means_B1_slow = PS_means_B1_slow * 1e4;
PS_means_B1_fast_SXL = PS_means_B1_fast_SXL * 1e4;
PS_means_B1_slow_SXL = PS_means_B1_slow_SXL * 1e4;
PS_means_B1_fast_corrected = PS_means_B1_fast_corrected* 1e4;
PS_means_B1_slow_corrected = PS_means_B1_slow_corrected* 1e4;
PS_means_B1_fast_SXL_corrected = PS_means_B1_fast_SXL_corrected* 1e4;
PS_means_B1_slow_SXL_corrected = PS_means_B1_slow_SXL_corrected* 1e4;
PS_devs_B1_fast = PS_devs_B1_fast * 1e4;
PS_devs_B1_slow = PS_devs_B1_slow * 1e4;
PS_devs_B1_fast_SXL = PS_devs_B1_fast_SXL * 1e4;
PS_devs_B1_slow_SXL = PS_devs_B1_slow_SXL * 1e4;
PS_devs_B1_fast_corrected = PS_devs_B1_fast_corrected * 1e4;
PS_devs_B1_slow_corrected = PS_devs_B1_slow_corrected * 1e4;
PS_devs_B1_fast_SXL_corrected = PS_devs_B1_fast_SXL_corrected * 1e4;
PS_devs_B1_slow_SXL_corrected = PS_devs_B1_slow_SXL_corrected * 1e4;

vP_range = vP_range * 1e2;
vP_means_B1_fast = vP_means_B1_fast * 1e2;
vP_means_B1_slow = vP_means_B1_slow * 1e2;
vP_means_B1_fast_SXL = vP_means_B1_fast_SXL * 1e2;
vP_means_B1_slow_SXL = vP_means_B1_slow_SXL * 1e2;
vP_means_B1_fast_corrected = vP_means_B1_fast_corrected* 1e2;
vP_means_B1_slow_corrected = vP_means_B1_slow_corrected* 1e2;
vP_means_B1_fast_SXL_corrected = vP_means_B1_fast_SXL_corrected* 1e2;
vP_means_B1_slow_SXL_corrected = vP_means_B1_slow_SXL_corrected* 1e2;
vP_devs_B1_fast = vP_devs_B1_fast * 1e2;
vP_devs_B1_slow = vP_devs_B1_slow * 1e2;
vP_devs_B1_fast_SXL = vP_devs_B1_fast_SXL * 1e2;
vP_devs_B1_slow_SXL = vP_devs_B1_slow_SXL * 1e2;
vP_devs_B1_fast_corrected = vP_devs_B1_fast_corrected * 1e2;
vP_devs_B1_slow_corrected = vP_devs_B1_slow_corrected * 1e2;
vP_devs_B1_fast_SXL_corrected = vP_devs_B1_fast_SXL_corrected * 1e2;
vP_devs_B1_slow_SXL_corrected = vP_devs_B1_slow_SXL_corrected * 1e2;

 save('PS_means_B1','PS_means_B1_fast','PS_means_B1_slow','PS_means_B1_fast_corrected','PS_means_B1_slow_corrected'...
     ,'PS_means_B1_fast_SXL','PS_means_B1_slow_SXL','PS_means_B1_fast_SXL_corrected','PS_means_B1_slow_SXL_corrected')
 save('PS_devs_B1','PS_devs_B1_fast','PS_devs_B1_slow','PS_devs_B1_fast_corrected','PS_devs_B1_slow_corrected'...
     ,'PS_devs_B1_fast_SXL','PS_devs_B1_slow_SXL','PS_devs_B1_fast_SXL_corrected','PS_devs_B1_slow_SXL_corrected')
 save('vP_means_B1','vP_means_B1_fast','vP_means_B1_slow','vP_means_B1_fast_corrected','vP_means_B1_slow_corrected'...
     ,'vP_means_B1_fast_SXL','vP_means_B1_slow_SXL','vP_means_B1_fast_SXL_corrected','vP_means_B1_slow_SXL_corrected')
 save('vP_devs_B1','vP_devs_B1_fast','vP_devs_B1_slow','vP_devs_B1_fast_corrected','vP_devs_B1_slow_corrected'...
     ,'vP_devs_B1_fast_SXL','vP_devs_B1_slow_SXL','vP_devs_B1_fast_SXL_corrected','vP_devs_B1_slow_SXL_corrected')
 
 % calculate average of PS or Vp error for a measure of
% systematic bias over the range of tested values
mean_PS_bias_fast_SXL_KFA112 = mean2(PS_means_B1_fast_SXL(:,2) - repmat(PS_range,1,1))
mean_PS_bias_slow_SXL_KFA112 = mean2(PS_means_B1_slow_SXL(:,2) - repmat(PS_range,1,1))
mean_PS_bias_fast_SXL_KFAdifferent = mean2(PS_means_B1_fast_SXL(:,3) - repmat(PS_range,1,1))
mean_PS_bias_slow_SXL_KFAdifferent = mean2(PS_means_B1_slow_SXL(:,3) - repmat(PS_range,1,1))

%% Calculate sensitivities
PS_sensitivities = NaN(1,8);
vP_sensitivities = NaN(1,8);
for i = 1:10
    PS_s_fast(i) = max(PS_means_B1_fast(i,:)) - min(PS_means_B1_fast(i,:));
    PS_s_fast_corrected(i) = max(PS_means_B1_fast_corrected(i,:)) - min(PS_means_B1_fast_corrected(i,:));
    PS_s_fast_SXL(i) = max(PS_means_B1_fast_SXL(i,:)) - min(PS_means_B1_fast_SXL(i,:));
    PS_s_fast_SXL_corrected(i) = max(PS_means_B1_fast_SXL_corrected(i,:)) - min(PS_means_B1_fast_SXL_corrected(i,:));
    PS_s_slow(i) = max(PS_means_B1_slow(i,:)) - min(PS_means_B1_slow(i,:));
    PS_s_slow_corrected(i) = max(PS_means_B1_slow_corrected(i,:)) - min(PS_means_B1_slow_corrected(i,:));
    PS_s_slow_SXL(i) = max(PS_means_B1_slow_SXL(i,:)) - min(PS_means_B1_slow_SXL(i,:));
    PS_s_slow_SXL_corrected(i) = max(PS_means_B1_slow_SXL_corrected(i,:)) - min(PS_means_B1_slow_SXL_corrected(i,:));
    
    vP_s_fast(i) = max(vP_means_B1_fast(i,:)) - min(vP_means_B1_fast(i,:));
    vP_s_fast_corrected(i) = max(vP_means_B1_fast_corrected(i,:)) - min(vP_means_B1_fast_corrected(i,:));
    vP_s_fast_SXL(i) = max(vP_means_B1_fast_SXL(i,:)) - min(vP_means_B1_fast_SXL(i,:));
    vP_s_fast_SXL_corrected(i) = max(vP_means_B1_fast_SXL_corrected(i,:)) - min(vP_means_B1_fast_SXL_corrected(i,:));
    vP_s_slow(i) = max(vP_means_B1_slow(i,:)) - min(vP_means_B1_slow(i,:));
    vP_s_slow_corrected(i) = max(vP_means_B1_slow_corrected(i,:)) - min(vP_means_B1_slow_corrected(i,:));
    vP_s_slow_SXL(i) = max(vP_means_B1_slow_SXL(i,:)) - min(vP_means_B1_slow_SXL(i,:));
    vP_s_slow_SXL_corrected(i) = max(vP_means_B1_slow_SXL_corrected(i,:)) - min(vP_means_B1_slow_SXL_corrected(i,:));
end    
    PS_sensitivities(1) = mean(PS_s_fast);
    PS_sensitivities(2) = mean(PS_s_slow);
    PS_sensitivities(3) = mean(PS_s_fast_SXL);
    PS_sensitivities(4) = mean(PS_s_slow_SXL);
    PS_sensitivities(5) = mean(PS_s_fast_corrected);
    PS_sensitivities(6) = mean(PS_s_slow_corrected);
    PS_sensitivities(7) = mean(PS_s_fast_SXL_corrected);
    PS_sensitivities(8) = mean(PS_s_slow_SXL_corrected);
    
    vP_sensitivities(1) = mean(vP_s_fast);
    vP_sensitivities(2) = mean(vP_s_slow);
    vP_sensitivities(3) = mean(vP_s_fast_SXL);
    vP_sensitivities(4) = mean(vP_s_slow_SXL);
    vP_sensitivities(5) = mean(vP_s_fast_corrected);
    vP_sensitivities(6) = mean(vP_s_slow_corrected);
    vP_sensitivities(7) = mean(vP_s_fast_SXL_corrected);
    vP_sensitivities(8) = mean(vP_s_slow_SXL_corrected);
    
%% Figure of PS values
Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];
%Colour4 = [0.4940 0.1840 0.5560 0.5];

figure()
h1=subplot(2,4,1) % fast, FXL, no correction

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_fast(:,1) - PS_range, 1*PS_devs_B1_fast(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_fast(:,2) - PS_range, 1*PS_devs_B1_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_fast(:,3) - PS_range, 1*PS_devs_B1_fast(:,3),'LineWidth',1.1,'Color',Colour3);
ylabel({'{\bfNo {\itB_1^+} correction}'},'FontSize',8);
%xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)]);
ylim([-2 2]);
title('Bolus (FXL fitting)','FontSize',8);

%annotation(figure(1),'textbox',[0.138 0.859 0.3 0.045],'String',['{\its} = ' num2str(round(PS_sensitivities(1),2)) ' x 10^{-4} min^{-1}'],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,2) % slow, FXL, no correction

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_slow(:,1) - PS_range, 1*PS_devs_B1_slow(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_slow(:,2) - PS_range, 1*PS_devs_B1_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_slow(:,3) - PS_range, 1*PS_devs_B1_slow(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-2 2]);
title('Slow inj (FXL fitting)','FontSize',8);
%xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);

%annotation(figure(1),'textbox',[0.345 0.584 0.3 0.045],'String',['{\its} = 1.10 x 10^{-4} min^{-1}'],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,3) % fast, NXL, no correction

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_fast_SXL(:,1) - PS_range, 1*PS_devs_B1_fast_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_fast_SXL(:,2) - PS_range, 1*PS_devs_B1_fast_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_fast_SXL(:,3) - PS_range, 1*PS_devs_B1_fast_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-2 2]);
title('Bolus (NXL fitting)','FontSize',8);
%xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);

%annotation(figure(1),'textbox',[0.552 0.859 0.3 0.045],'String',['{\its} = ' num2str(round(PS_sensitivities(3),2)) ' x 10^{-4} min^{-1}'],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,4) % slow, NXL, no correction

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_slow_SXL(:,1) - PS_range, 1*PS_devs_B1_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_slow_SXL(:,2) - PS_range, 1*PS_devs_B1_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_slow_SXL(:,3) - PS_range, 1*PS_devs_B1_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-2 2]);
title('Slow inj (NXL fitting)','FontSize',8);
%xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);

%annotation(figure(1),'textbox',[0.759 0.859 0.3 0.045],'String',['{\its} = 0.90 x 10^{-4} min^{-1}'],'LineStyle','none','FitBoxToText','off','FontSize',8);

h2=subplot(2,4,5) % fast, FXL, correction

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_fast_corrected(:,1) - PS_range, 1*PS_devs_B1_fast_corrected(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_fast_corrected(:,2) - PS_range, 1*PS_devs_B1_fast_corrected(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_fast_corrected(:,3) - PS_range, 1*PS_devs_B1_fast_corrected(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-0.8 2]);
ylabel({'{\bf{\itB_1^+} corrected}'},'FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
legend({'1,1', '1.12,1.12','0.95,1.12'},'Position',[0.172 0.360 0.081 0.0207],'FontSize',8)
legend('boxoff')

%annotation(figure(1),'textbox',[0.138 0.110 0.3 0.045],'String',['{\its} = ' num2str(round(PS_sensitivities(5),2)) ' x 10^{-4} min^{-1}'],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,6) % slow, FXL, correction

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_slow_corrected(:,1) - PS_range, 1*PS_devs_B1_slow_corrected(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_slow_corrected(:,2) - PS_range, 1*PS_devs_B1_slow_corrected(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_slow_corrected(:,3) - PS_range, 1*PS_devs_B1_slow_corrected(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-0.8 2]);
%title('Slow injection','FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);

%annotation(figure(1),'textbox',[0.345 0.110 0.3 0.045],'String',['{\its} = ' num2str(round(PS_sensitivities(6),2)) ' x 10^{-4} min^{-1}'],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,7) % fast, SXL, correction

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_fast_SXL_corrected(:,1) - PS_range, 1*PS_devs_B1_fast_SXL_corrected(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_fast_SXL_corrected(:,2) - PS_range, 1*PS_devs_B1_fast_SXL_corrected(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_fast_SXL_corrected(:,3) - PS_range, 1*PS_devs_B1_fast_SXL_corrected(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-0.8 2]);
%title('Slow injection','FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);

%annotation(figure(1),'textbox',[0.552 0.385 0.3 0.045],'String',['{\its} = ' num2str(round(PS_sensitivities(7),2)) ' x 10^{-4} min^{-1}'],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,8) % slow, SXL, correction

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_slow_SXL_corrected(:,1) - PS_range, 1*PS_devs_B1_slow_SXL_corrected(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_slow_SXL_corrected(:,2) - PS_range, 1*PS_devs_B1_slow_SXL_corrected(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_slow_SXL_corrected(:,3) - PS_range, 1*PS_devs_B1_slow_SXL_corrected(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-0.8 2]);
%title('Slow injection','FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);

%annotation(figure(1),'textbox',[0.759 0.385 0.3 0.045],'String',['{\its} = ' num2str(round(PS_sensitivities(8),2)) ' x 10^{-4} min^{-1}'],'LineStyle','none','FitBoxToText','off','FontSize',8);

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

annotation(figure(1),'textbox',[0.155 0.405 0.2 0.045],'String','{\itK_{FA,b}} , {\itK_{FA,t}}','LineStyle','none','FitBoxToText','off','FontSize',9);

set(gcf, 'units', 'centimeters','PaperPosition', [0 0 20 16]);    % can be bigger than screen
print(gcf, 'Figure_5.png', '-dpng','-r800');
print(gcf, 'Figure_5.tif', '-dtiff','-r800');

%% Figure of vp values
figure(2)

h3=subplot(2,4,1) % fast, FXL, no correction

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_fast(:,1) - vP_range, 1*vP_devs_B1_fast(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_B1_fast(:,2) - vP_range, 1*vP_devs_B1_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_B1_fast(:,3) - vP_range, 1*vP_devs_B1_fast(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
%xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
ylabel({'{\bfNo {\itB_1^+} correction}'},'FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-1 0.2]);
title('Bolus (FXL fitting)','FontSize',8);

%annotation(figure(2),'textbox',[0.168 0.584 0.3 0.045],'String',['{\its} = ' num2str(round(vP_sensitivities(1),2))],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,2) % slow, FXL, no correction

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_slow(:,1) - vP_range, 1*vP_devs_B1_slow(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_B1_slow(:,2) - vP_range, 1*vP_devs_B1_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_B1_slow(:,3) - vP_range, 1*vP_devs_B1_slow(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
%xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-1 0.2]);
title('Slow inj (FXL fitting)','FontSize',8);

%annotation(figure(2),'textbox',[0.375 0.584 0.3 0.045],'String',['{\its} = ' num2str(round(vP_sensitivities(2),2))],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,3) % fast, NXL, no correction

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_fast_SXL(:,1) - vP_range, 1*vP_devs_B1_fast_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_B1_fast_SXL(:,2) - vP_range, 1*vP_devs_B1_fast_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_B1_fast_SXL(:,3) - vP_range, 1*vP_devs_B1_fast_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
%xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-1 0.2]);
title('Bolus (NXL fitting)','FontSize',8);

%annotation(figure(2),'textbox',[0.582 0.584 0.3 0.045],'String',['{\its} = ' num2str(round(vP_sensitivities(3),2))],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,4) % slow, NXL, no correction

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_slow_SXL(:,1) - vP_range, 1*vP_devs_B1_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_B1_slow_SXL(:,2) - vP_range, 1*vP_devs_B1_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_B1_slow_SXL(:,3) - vP_range, 1*vP_devs_B1_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
%xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-1 0.2]);
title('Slow inj (NXL fitting)','FontSize',8);

%annotation(figure(2),'textbox',[0.789 0.584 0.3 0.045],'String',['{\its} = 0.30'],'LineStyle','none','FitBoxToText','off','FontSize',8);

h4=subplot(2,4,5) % fast, FXL, correction

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_fast_corrected(:,1) - vP_range, 1*vP_devs_B1_fast_corrected(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_B1_fast_corrected(:,2) - vP_range, 1*vP_devs_B1_fast_corrected(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_B1_fast_corrected(:,3) - vP_range, 1*vP_devs_B1_fast_corrected(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
ylabel({'{\bf{\itB_1^+} corrected}'},'FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.2]);

%annotation(figure(2),'textbox',[0.168 0.385 0.3 0.045],'String',['{\its} = ' num2str(round(vP_sensitivities(5),2))],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,6) % slow, FXL, correction

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_slow_corrected(:,1) - vP_range, 1*vP_devs_B1_slow_corrected(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_B1_slow_corrected(:,2) - vP_range, 1*vP_devs_B1_slow_corrected(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_B1_slow_corrected(:,3) - vP_range, 1*vP_devs_B1_slow_corrected(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.2]);

%annotation(figure(2),'textbox',[0.375 0.385 0.3 0.045],'String',['{\its} = ' num2str(round(vP_sensitivities(6),2))],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,7) % fast, NXL, correction

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_fast_SXL_corrected(:,1) - vP_range, 1*vP_devs_B1_fast_SXL_corrected(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_B1_fast_SXL_corrected(:,2) - vP_range, 1*vP_devs_B1_fast_SXL_corrected(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_B1_fast_SXL_corrected(:,3) - vP_range, 1*vP_devs_B1_fast_SXL_corrected(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.2]);
legend({'1,1', '1.12,1.12', '0.95,1.12'},'Location','best','FontSize',8)
legend('boxoff')

%annotation(figure(2),'textbox',[0.582 0.385 0.3 0.045],'String',['{\its} = ' num2str(round(vP_sensitivities(7),2))],'LineStyle','none','FitBoxToText','off','FontSize',8);

subplot(2,4,8) % slow, NXL, correction

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_slow_SXL_corrected(:,1) - vP_range, 1*vP_devs_B1_slow_SXL_corrected(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.013, vP_means_B1_slow_SXL_corrected(:,2) - vP_range, 1*vP_devs_B1_slow_SXL_corrected(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.026, vP_means_B1_slow_SXL_corrected(:,3) - vP_range, 1*vP_devs_B1_slow_SXL_corrected(:,3),'LineWidth',1.1,'Color',Colour3); hold on;
xlabel('True {\itv_p} (x10^{-2})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.026]);
ylim([-0.5 0.2]);

%annotation(figure(2),'textbox',[0.789 0.110 0.3 0.045],'String',['{\its} = ' num2str(round(vP_sensitivities(8),2))],'LineStyle','none','FitBoxToText','off','FontSize',8);

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

annotation(figure(2),'textbox',[0.563 0.184 0.2 0.045],'String','{\itK_{FA,b}} , {\itK_{FA,t}}','LineStyle','none','FitBoxToText','off','FontSize',9);

set(gcf, 'units', 'centimeters','PaperPosition', [0 0 20 16]);    % can be bigger than screen
set(gcf,'PaperPositionMode','manual');
print(gcf, 'Supp_Figure_vP_B1.png', '-dpng','-r800');