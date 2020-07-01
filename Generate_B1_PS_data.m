% Generates the B1 inhomogeneity PS simulation data as shown in Manning et al.
% (2020) Slow injection paper
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;

addpath('DCE_Simulation_Functions');
    
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.water_exch_model = '2S1XA';
SimParam.SXLfit = 0; % fit enhancements according to SXL method

PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;

N_PS = size(PS_range,1); %range sizes to test
%% Generate B1 inhomogeneity sims
%% B1 inhomogeneity figures (fast injection, no exclude)

 for i_PS = 1:N_PS % Accurate T1 acquisition
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_highkaccurate_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,1) = mean(PS_fit_highkaccurate_fast,1)'; % mean for each PS for k = 1.3, double correction
PS_devs_VFA_fast(:,1) = std(PS_fit_highkaccurate_fast,0,1)'; % standard deviation at k = 1.3, double correction

T1acqParam.T1_acq_method = 'VFA';
DCESeqParam.FA_error = 1.3;
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
T1acqParam.T1_SNR = 319; % to achieve comparable T1 error in WM (Lee, Callaghan, et al 2018)
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS % VFA T1 acquisition, k = 1.3
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_highkVFA_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,2) = mean(PS_fit_highkVFA_fast,1)'; % mean for each PS for k = 1.3, uncorrected
PS_devs_VFA_fast(:,2) = std(PS_fit_highkVFA_fast,0,1)'; % standard deviation at k = 1.3, uncorrected

T1acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 0.7;
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
T1acqParam.T1_SNR = 319; % to achieve comparable T1 error in WM (Lee, Callaghan, et al 2018)
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

 for i_PS = 1:N_PS % VFA T1 acquisition, k = 0.7
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_lowkVFA_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_fast(:,3) = mean(PS_fit_lowkVFA_fast,1)'; % mean for each PS for k = 0.7, uncorrected
PS_devs_VFA_fast(:,3) = std(PS_fit_lowkVFA_fast,0,1)'; % standard deviation at k = 0.7, uncorrected

%% Generate variable flow PS and vP - fast injection (exclude)
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.water_exch_model = '2S1XA';
SimParam.SXLfit = 0; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;
  
 for i_PS = 1:N_PS % Accurate T1 acquisition
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_highkaccurate_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_exclude(:,1) = mean(PS_fit_highkaccurate_exclude,1)'; % mean for each PS for k = 1.3, double correction
PS_devs_VFA_exclude(:,1) = std(PS_fit_highkaccurate_exclude,0,1)'; % standard deviation at k = 1.3, double correction

% Simulate T1 acquisiton
T1acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.3;
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
T1acqParam.T1_SNR = 319; % to achieve comparable T1 error in WM (Lee, Callaghan, et al 2018)
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS % VFA T1 acquisition, k = 1.3
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_highkVFA_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_exclude(:,2) = mean(PS_fit_highkVFA_exclude,1)'; % mean for each PS for k = 1.3
PS_devs_VFA_exclude(:,2) = std(PS_fit_highkVFA_exclude,0,1)'; % standard deviation for k = 1.3

DCESeqParam.FA_error = 0.7;
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
T1acqParam.T1_SNR = inf;
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
T1acqParam.T1_SNR = 319; % to achieve comparable T1 error in WM (Lee, Callaghan, et al 2018)
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

for i_PS = 1:N_PS % VFA T1 acquisition, k = 0.7
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_lowkVFA_exclude(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
end
 
PS_means_VFA_exclude(:,3) = mean(PS_fit_lowkVFA_exclude,1)'; % mean for each PS for k = 0.7
PS_devs_VFA_exclude(:,3) = std(PS_fit_lowkVFA_exclude,0,1)'; % standard deviation at k = 0.7

%% Generate B1 inhomogeneity figures - slow injection
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.water_exch_model = '2S1XA';
SimParam.SXLfit = 0; % fit enhancements according to SXL method

 SimParam.t_start_s = 0;
 SimParam.InjectionRate = 'slow';
 SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points

 for i_PS = 1:N_PS % Accurate T1 measurement
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_highkaccurate_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_slow(:,1) = mean(PS_fit_highkaccurate_slow,1)'; % mean for each PS for k = 1.3, double correction
PS_devs_VFA_slow(:,1) = std(PS_fit_highkaccurate_slow,0,1)'; % standard deviation at k = 1.3, double correction

T1acqParam.T1_acq_method = 'VFA'; % Accurate T1 measurement
DCESeqParam.FA_error = 1.3;
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
T1acqParam.T1_SNR = 319; % to achieve comparable T1 error in WM (Lee, Callaghan, et al 2018)
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_nom_deg; % true flip angle

 for i_PS = 1:N_PS  % VFA T1 acquisition, k = 1.3
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_highkaccurate_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
 end
 
PS_means_VFA_slow(:,2) = mean(PS_fit_highkaccurate_slow,1)'; % mean for each PS for k = 0.7, double correction
PS_devs_VFA_slow(:,2) = std(PS_fit_highkaccurate_slow,0,1)'; % standard deviation at k = 0.7, double correction

% Simulate T1 acquisiton
DCESeqParam.FA_error = 0.7; 
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
T1acqParam.T1_SNR = 319; % to achieve comparable T1 error in WM (Lee, Callaghan, et al 2018)
[PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle

   for i_PS = 1:N_PS % VFA T1 acquisition, k = 0.7
    PhysParam.vP = vP_fixed(1);
    PhysParam.PS_perMin = PS_range(i_PS);
    [temp, PS_fit_lowkVFA_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
  end
  
PS_means_VFA_slow(:,3) = mean(PS_fit_lowkVFA_slow,1)'; % add mean for each PS for k = 1.3, uncorrected
PS_devs_VFA_slow(:,3) = std(PS_fit_lowkVFA_slow,0,1)'; % add standard deviation for k = 1.3, uncorrected

%% Graphs and scales
PS_range = PS_range * 1e4;
PS_means_VFA_fast = PS_means_VFA_fast * 1e4;
PS_means_VFA_exclude = PS_means_VFA_exclude * 1e4;
PS_means_VFA_slow = PS_means_VFA_slow * 1e4;
PS_devs_VFA_fast = PS_devs_VFA_fast * 1e4;
PS_devs_VFA_exclude = PS_devs_VFA_exclude * 1e4;
PS_devs_VFA_slow = PS_devs_VFA_slow * 1e4;

 save('PS_means_VFA','PS_means_VFA_fast','PS_means_VFA_exclude','PS_means_VFA_slow')
 save('PS_devs_VFA','PS_devs_VFA_fast','PS_devs_VFA_exclude','PS_devs_VFA_slow')

Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];
%Colour4 = [0.4940 0.1840 0.5560 0.5];

subplot(1,3,1)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_fast(:,1) - PS_range, 1*PS_devs_VFA_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_fast(:,2) - PS_range, 1*PS_devs_VFA_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_fast(:,3) - PS_range, 1*PS_devs_VFA_fast(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
%errorbar(PS_range + 0.09, PS_means_VFA_fast(:,4) - PS_range, 1*PS_devs_VFA_fast(:,4),'LineWidth',1.3,'Color',Colour4); hold on;
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-4 4]);
legend({'Accurate T_1_0 and DCE','T_1_0 and DCE uncorrected (k=1.3)','T_1_0 and DCE uncorrected (k=0.7)'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 9;

subplot(1,3,2)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_exclude(:,1) - PS_range, 1*PS_devs_VFA_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_exclude(:,2) - PS_range, 1*PS_devs_VFA_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_exclude(:,3) - PS_range, 1*PS_devs_VFA_exclude(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
%errorbar(PS_range + 0.09, PS_means_VFA_exclude(:,4) - PS_range, 1*PS_devs_VFA_exclude(:,4),'LineWidth',1.3,'Color',Colour4);
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(1,3,3)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_slow(:,1) - PS_range, 1*PS_devs_VFA_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_slow(:,2) - PS_range, 1*PS_devs_VFA_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_slow(:,3) - PS_range, 1*PS_devs_VFA_slow(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
%errorbar(PS_range + 0.09, PS_means_VFA_slow(:,4) - PS_range, 1*PS_devs_VFA_slow(:,4),'LineWidth',1.3,'Color',Colour4); hold on;
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

set(gcf, 'units', 'centimeters','Position', [5 5 28 10]);