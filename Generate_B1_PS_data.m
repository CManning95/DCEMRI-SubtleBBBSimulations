% Generates the B1 inhomogeneity PS simulation data as shown in Manning et al.
% (2020) Slow injection paper
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;

addpath('DCE_Simulation_Functions');
    
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.SXLfit = 1; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;

PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;
vP_range = linspace(SimParam.min_vP,SimParam.max_vP,10)';
PS_fixed = PhysParam.PS_fixed;
k_values = [1 0.7 1.3]; 

N_PS = size(PS_range,1); %range sizes to test
N_vP = size(vP_range,1); %range sizes to test

%% Generate B1 inhomogeneity sims
%% B1 inhomogeneity figures (fast injection)
for k = 1:size(k_values,2)
    DCESeqParam.FA_error = k_values(k);
    T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
    DCESeqParam.FA_true_deg = DCESeqParam.FA_error*DCESeqParam.FA_meas_deg; % true flip angle
    [PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
    DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
    
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    for i_vP = 1:N_vP
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
%% B1 inhomogeneity figures (fast injection - B1 corrected)
for k = 1:size(k_values,2)
    T1acqParam.T1_acq_method = 'HIFI';
    T1acqParam.isFit = [1 1 1 1 1]; % which acquisitions to fit
    T1acqParam.TR_s = [0.0054 0.0054 0.0054 0.0054 0.0054]; % repetition times for T1 acqusition
    T1acqParam.FA_nom_rads = [5 5 2 5 12] *2*(pi/360); % nominal flip angles for T1 acquisition
    T1acqParam.isIR = [1 1 0 0 0]; % indicates which are IR-SPGR
    T1acqParam.TI_s = [0.168 1.068 NaN NaN NaN]; % Inversion times for T1 acquisition (for HIFI)
    T1acqParam.PECentre = [0.5 0.5 NaN NaN NaN]; % indicates time of centre of k-space
    T1acqParam.NReadout = [160 160 160 160 160]; % number of readout pulses (Siemens - number of slices)
    
    DCESeqParam.FA_error = k_values(k);
    T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
    DCESeqParam.FA_true_deg = DCESeqParam.FA_error * DCESeqParam.FA_nom_deg; % true flip angle
    [PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
    DCESeqParam.FA_meas_deg = T1acqParam.FA_error_meas * DCESeqParam.FA_nom_deg;
    
    for i_PS = 1:N_PS % Accurate T1 acquisition
        for n = 1:SimParam.N_repetitions
            [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
        end
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_fast_corrected(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    for i_vP = 1:N_vP
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

%% Generate variable flow PS and vP - slow injection
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.SXLfit = 1; % fit enhancements according to SXL method
SimParam.NIgnore = max(SimParam.baselineScans) + 3;
 SimParam.t_start_s = 0;
 SimParam.InjectionRate = 'slow';
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
 
 for k = 1:size(k_values,2)
     DCESeqParam.FA_error = k_values(k);
     T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
     DCESeqParam.FA_true_deg = DCESeqParam.FA_error * DCESeqParam.FA_nom_deg; % true flip angle
     [PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
     DCESeqParam.FA_meas_deg = DCESeqParam.FA_nom_deg;
     
     for i_PS = 1:N_PS % Accurate T1 measurement
         for n = 1:SimParam.N_repetitions
             [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
         end
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
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

%% Generate B1 inhomogeneity figures - (slow injection - B1 corrected)
for k = 1:size(k_values,2)
    T1acqParam.T1_acq_method = 'HIFI';
    T1acqParam.isFit = [1 1 1 1 1]; % which acquisitions to fit
    T1acqParam.TR_s = [0.0054 0.0054 0.0054 0.0054 0.0054]; % repetition times for T1 acqusition
    T1acqParam.FA_nom_rads = [5 5 2 5 12] *2*(pi/360); % nominal flip angles for T1 acquisition
    T1acqParam.isIR = [1 1 0 0 0]; % indicates which are IR-SPGR
    T1acqParam.TI_s = [0.168 1.068 NaN NaN NaN]; % Inversion times for T1 acquisition (for HIFI)
    T1acqParam.PECentre = [0.5 0.5 NaN NaN NaN]; % indicates time of centre of k-space
    T1acqParam.NReadout = [160 160 160 160 160]; % number of readout pulses (Siemens - number of slices)
    
    DCESeqParam.FA_error = k_values(k);
    T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
    DCESeqParam.FA_true_deg = DCESeqParam.FA_error * DCESeqParam.FA_nom_deg; % true flip angle
    [PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
    [PhysParam.T1_tissue_meas_s,temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
    DCESeqParam.FA_meas_deg = T1acqParam.FA_error_meas*DCESeqParam.FA_nom_deg;
    
    for i_PS = 1:N_PS % Accurate T1 acquisition
         for n = 1:SimParam.N_repetitions
             [T1_tissue_meas_s(n,i_PS),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
         end
        PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
        PhysParam.vP = vP_fixed(1);
        PhysParam.PS_perMin = PS_range(i_PS);
        [temp, PS_fit_slow_corrected(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);    
    end
    for i_vP = 1:N_vP
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_vP);
         PhysParam.PS = PS_fixed(1);
         PhysParam.vP = vP_range(i_vP);
         [vP_fit_slow_corrected(:,i_vP),temp] = master_single_sim(PhysParam,DCESeqParam,SimParam);
    end
    
    PS_means_B1_slow_corrected(:,k) = mean(PS_fit_slow_corrected,1)'; % mean for each PS
    PS_devs_B1_slow_corrected(:,k) = std(PS_fit_slow_corrected,0,1)'; % standard deviation for each PS
    
    vP_means_B1_slow_corrected(:,k) = mean(vP_fit_slow_corrected,1)'; % mean for each PS
    vP_devs_B1_slow_corrected(:,k) = std(vP_fit_slow_corrected,0,1)'; % standard deviation for each PS
end

%% Graphs and scales
PS_range = PS_range * 1e4;
PS_means_B1_fast = PS_means_B1_fast * 1e4;
PS_means_B1_slow = PS_means_B1_slow * 1e4;
PS_means_B1_fast_corrected = PS_means_B1_fast_corrected* 1e4;
PS_means_B1_slow_corrected = PS_means_B1_slow_corrected* 1e4;
PS_devs_B1_fast = PS_devs_B1_fast * 1e4;
PS_devs_B1_slow = PS_devs_B1_slow * 1e4;
PS_devs_B1_fast_corrected = PS_devs_B1_fast_corrected * 1e4;
PS_devs_B1_slow_corrected = PS_devs_B1_slow_corrected * 1e4;

vP_range = vP_range * 1e3;
vP_means_B1_fast = vP_means_B1_fast * 1e3;
vP_means_B1_slow = vP_means_B1_slow * 1e3;
vP_means_B1_fast_corrected = vP_means_B1_fast_corrected* 1e3;
vP_means_B1_slow_corrected = vP_means_B1_slow_corrected* 1e3;
vP_devs_B1_fast = vP_devs_B1_fast * 1e3;
vP_devs_B1_slow = vP_devs_B1_slow * 1e3;
vP_devs_B1_fast_corrected = vP_devs_B1_fast_corrected * 1e3;
vP_devs_B1_slow_corrected = vP_devs_B1_slow_corrected * 1e3;

%  save('PS_means_B1','PS_means_B1_fast','PS_means_B1_slow','PS_means_B1_fast_corrected','PS_means_B1_slow_corrected')
%  save('PS_devs_B1','PS_devs_B1_fast','PS_devs_B1_slow','PS_devs_B1_fast_corrected','PS_devs_B1_slow_corrected')
%  save('vP_means_B1','vP_means_B1_fast','vP_means_B1_slow','vP_means_B1_fast_corrected','vP_means_B1_slow_corrected')
%  save('vP_devs_B1','vP_devs_B1_fast','vP_devs_B1_slow','vP_devs_B1_fast_corrected','vP_devs_B1_slow_corrected')
%%
Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];
%Colour4 = [0.4940 0.1840 0.5560 0.5];

figure()
subplot(2,4,1)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_fast(:,1) - PS_range, 1*PS_devs_B1_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_fast(:,2) - PS_range, 1*PS_devs_B1_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_fast(:,3) - PS_range, 1*PS_devs_B1_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-2 2]);
title(['Bolus injection']);
legend({'k = 1', 'k = 0.7', 'k = 1.3'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 9;

subplot(2,4,2)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_slow(:,1) - PS_range, 1*PS_devs_B1_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_slow(:,2) - PS_range, 1*PS_devs_B1_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_slow(:,3) - PS_range, 1*PS_devs_B1_slow(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-2 2]);
title(['Slow injection']);

ax = gca;
ax.FontSize = 9;

subplot(2,4,3)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_fast_corrected(:,1) - PS_range, 1*PS_devs_B1_fast_corrected(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_fast_corrected(:,2) - PS_range, 1*PS_devs_B1_fast_corrected(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_fast_corrected(:,3) - PS_range, 1*PS_devs_B1_fast_corrected(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-2 2]);
title(['Bolus (B1 corrected)']);

subplot(2,4,4)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_B1_slow_corrected(:,1) - PS_range, 1*PS_devs_B1_slow_corrected(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_B1_slow_corrected(:,2) - PS_range, 1*PS_devs_B1_slow_corrected(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_B1_slow_corrected(:,3) - PS_range, 1*PS_devs_B1_slow_corrected(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
xlim([0 max(PS_range)]);
ylim([-2 2]);
title(['Slow (B1 corrected)']);

subplot(2,4,5)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_fast(:,1) - vP_range, 1*vP_devs_B1_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_B1_fast(:,2) - vP_range, 1*vP_devs_B1_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_B1_fast(:,3) - vP_range, 1*vP_devs_B1_fast(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
ylabel('fitted v_p error (x10^{-3})');
xlabel(['True v_p (x10^{-3})']);
xlim([min(vP_range) max(vP_range)+0.26]);
ylim([-2 2]);

subplot(2,4,6)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_slow(:,1) - vP_range, 1*vP_devs_B1_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_B1_slow(:,2) - vP_range, 1*vP_devs_B1_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_B1_slow(:,3) - vP_range, 1*vP_devs_B1_slow(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
ylabel('fitted v_p error (x10^{-3})');
xlabel(['True v_p (x10^{-3})']);
xlim([min(vP_range) max(vP_range)+0.26]);
ylim([-2 2]);

subplot(2,4,7)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_fast_corrected(:,1) - vP_range, 1*vP_devs_B1_fast_corrected(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_B1_fast_corrected(:,2) - vP_range, 1*vP_devs_B1_fast_corrected(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_B1_fast_corrected(:,3) - vP_range, 1*vP_devs_B1_fast_corrected(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
ylabel('fitted v_p error (x10^{-3})');
xlabel(['True v_p (x10^{-3})']);
xlim([min(vP_range) max(vP_range)+0.26]);
ylim([-2 2]);

subplot(2,4,8)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True vP','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_B1_slow_corrected(:,1) - vP_range, 1*vP_devs_B1_slow_corrected(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_B1_slow_corrected(:,2) - vP_range, 1*vP_devs_B1_slow_corrected(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_B1_slow_corrected(:,3) - vP_range, 1*vP_devs_B1_slow_corrected(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
ylabel('fitted v_p error (x10^{-3})');
xlabel(['True v_p (x10^{-3})']);
xlim([min(vP_range) max(vP_range)+0.26]);
ylim([-2 2]);

set(gcf, 'units', 'centimeters','Position', [5 5 28 15]);

annotation(figure(1),'textbox',[0.084 0.935 0.05 0.045],'String','a. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.290 0.935 0.06 0.045],'String','a. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.494 0.935 0.06 0.045],'String','a. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.702 0.935 0.06 0.045],'String','a. iv','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.084 0.460 0.06 0.045],'String','b. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.290 0.460 0.06 0.045],'String','b. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.494 0.460 0.06 0.045],'String','b. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.702 0.460 0.06 0.045],'String','b. iv','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
set(gcf, 'units', 'centimeters','PaperPosition', [0 0 25 15]);    % can be bigger than screen 
print(gcf, 'B1_figure.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi