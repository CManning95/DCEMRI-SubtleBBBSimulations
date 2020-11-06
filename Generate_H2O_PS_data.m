
% Generates the water exchange PS simulation data as shown in Manning et al.
% (2020) Slow injection paper
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;

addpath('DCE_Simulation_Functions');
    
% Select default parameters
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
 
%%% Generate variable flow/injection delay sims
% ranges of PS and vP to test
PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;
vP_range = linspace(SimParam.min_vP,SimParam.max_vP,10)';
PS_fixed = PhysParam.PS_fixed;

N_PS = size(PS_range,1); %range sizes to test
N_vP = size(vP_range,1); %range sizes to test
kbe_ranges = [1.375 2.75 5.5];

%% Sim water exchange with Patlak fitting (fast injection, Patlak fit)
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
        
        PS_means_H2O_fast(:,i) = mean(PS_fit_2S1X_fast,1)'; % add mean for each PS for 2S1X
        PS_devs_H2O_fast(:,i) = std(PS_fit_2S1X_fast,0,1)'; % add standard deviation for 2S1X
        
        vP_means_H2O_fast(:,i) = mean(vP_fit_2S1X_fast,1)'; % add mean for each PS for 2S1X
        vP_devs_H2O_fast(:,i) = std(vP_fit_2S1X_fast,0,1)'; % add standard deviation for 2S1X
    end
    
    %% Sim water exchange with Patlak fitting (fast injection, SXL fit)
    SimParam.SXLfit = 1; % fit enhancements according to SXL method
    
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
        PS_means_H2O_SXL(:,i) = mean(PS_fit_2S1X_SXL,1)'; % add mean for each PS for 2S1X
        PS_devs_H2O_SXL(:,i) = std(PS_fit_2S1X_SXL,0,1)'; % add standard deviation for 2S1X
        
        vP_means_H2O_SXL(:,i) = mean(vP_fit_2S1X_SXL,1)'; % add mean for each PS for 2S1X
        vP_devs_H2O_SXL(:,i) = std(vP_fit_2S1X_SXL,0,1)'; % add standard deviation for 2S1X
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
        PS_means_H2O_slow(:,i) = mean(PS_fit_2S1X_slow,1)'; % add mean for each PS for 2S1X
        PS_devs_H2O_slow(:,i) = std(PS_fit_2S1X_slow,0,1)'; % add standard deviation for 2S1X
        
        vP_means_H2O_slow(:,i) = mean(vP_fit_2S1X_slow,1)'; % add mean for each PS for 2S1X
        vP_devs_H2O_slow(:,i) = std(vP_fit_2S1X_slow,0,1)'; % add standard deviation for 2S1X
        
    end
    
    %% Sim water exchange (slow injection, SXL fit)
    SimParam.SXLfit = 1;
    
    acqParam.T1_SNR = 318;
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
        PS_means_H2O_slow_SXL(:,i) = mean(PS_fit_2S1X_slow_SXL,1)'; % add mean for each PS for 2S1X
        PS_devs_H2O_slow_SXL(:,i) = std(PS_fit_2S1X_slow_SXL,0,1)'; % add standard deviation for 2S1X
        
        vP_means_H2O_slow_SXL(:,i) = mean(vP_fit_2S1X_slow,1)'; % add mean for each PS for 2S1X
        vP_devs_H2O_slow_SXL(:,i) = std(vP_fit_2S1X_slow,0,1)'; % add standard deviation for 2S1X
    end


    %% Save simulation data
    PS_range = PS_range * 1e4;
    PS_means_H2O_fast = PS_means_H2O_fast * 1e4;
    PS_means_H2O_SXL = PS_means_H2O_SXL * 1e4;
    PS_means_H2O_slow = PS_means_H2O_slow * 1e4;
    PS_means_H2O_slow_SXL = PS_means_H2O_slow_SXL * 1e4;
    PS_devs_H2O_fast = PS_devs_H2O_fast * 1e4;
    PS_devs_H2O_SXL = PS_devs_H2O_SXL * 1e4;
    PS_devs_H2O_slow = PS_devs_H2O_slow * 1e4;
    PS_devs_H2O_slow_SXL = PS_devs_H2O_slow_SXL * 1e4;
    save('PS_means_H2O','PS_means_H2O_fast','PS_means_H2O_SXL',...
        'PS_means_H2O_slow','PS_means_H2O_slow_SXL');
    save('PS_devs_H2O','PS_devs_H2O_fast','PS_devs_H2O_SXL',...
        'PS_devs_H2O_slow','PS_devs_H2O_slow_SXL');
    
    vP_range = vP_range * 1e3;
    vP_means_H2O_fast = vP_means_H2O_fast * 1e3;
    vP_means_H2O_SXL = vP_means_H2O_SXL * 1e3;
    vP_means_H2O_slow = vP_means_H2O_slow * 1e3;
    vP_means_H2O_slow_SXL = vP_means_H2O_slow_SXL * 1e3;
    vP_devs_H2O_fast = vP_devs_H2O_fast * 1e3;
    vP_devs_H2O_SXL = vP_devs_H2O_SXL * 1e3;
    vP_devs_H2O_slow = vP_devs_H2O_slow * 1e3;
    vP_devs_H2O_slow_SXL = vP_devs_H2O_slow_SXL * 1e3;
    save('vP_means_H2O','vP_means_H2O_fast','vP_means_H2O_SXL',...
        'vP_means_H2O_slow','vP_means_H2O_slow_SXL');
    save('vP_devs_H2O','vP_devs_H2O_fast','vP_devs_H2O_SXL',...
       'vP_devs_H2O_slow','vP_devs_H2O_slow_SXL');

Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];

%% Plot figures        
figure(1)

subplot(2,4,1)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_fast(:,1) - PS_range, 1*PS_devs_H2O_fast(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_H2O_fast(:,2) - PS_range, 1*PS_devs_H2O_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_fast(:,3) - PS_range, 1*PS_devs_H2O_fast(:,3),'LineWidth',1.1,'Color',Colour3);
ylabel('fitted {\itPS} error (x10^{-4} min^{-1} )','FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
title('Bolus (FXL fitting)','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-5 5]);
legend({'k_{be} = 1.375 s^{-1}','k_{be} = 2.75 s^{-1}','k_{be} = 5.5 s^{-1}'},'Location','best','FontSize',6)
legend('boxoff')

subplot(2,4,3)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_SXL(:,1) - PS_range, 1*PS_devs_H2O_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_H2O_SXL(:,2) - PS_range, 1*PS_devs_H2O_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_SXL(:,3) - PS_range, 1*PS_devs_H2O_SXL(:,3),'LineWidth',1.1,'Color',Colour3);
title('Bolus (SXL fitting)','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2 2]);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);

subplot(2,4,2)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_slow(:,1) - PS_range, 1*PS_devs_H2O_slow(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_H2O_slow(:,2) - PS_range, 1*PS_devs_H2O_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_slow(:,3) - PS_range, 1*PS_devs_H2O_slow(:,3),'LineWidth',1.1,'Color',Colour3);
title('Slow (FXL fitting)','FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2 2]);

subplot(2,4,4)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_slow_SXL(:,1) - PS_range, 1*PS_devs_H2O_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_H2O_slow_SXL(:,2) - PS_range, 1*PS_devs_H2O_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_H2O_slow_SXL(:,3) - PS_range, 1*PS_devs_H2O_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3);
title('Slow (SXL fitting)','FontSize',8);
xlabel('True {\itPS} (x10^{-4} min^{-1} )','FontSize',8);
xlim([0 max(PS_range)+0.12]);
ylim([-2 2]);

subplot(2,4,5)

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_fast(:,1) - vP_range, 1*vP_devs_H2O_fast(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_H2O_fast(:,2) - vP_range, 1*vP_devs_H2O_fast(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_H2O_fast(:,3) - vP_range, 1*vP_devs_H2O_fast(:,3),'LineWidth',1.1,'Color',Colour3);
xlabel('True {\itv_p} (x10^{-3})','FontSize',8);
ylabel('fitted {\itv_p} error (x10^{-3})','FontSize',8);
xlim([min(vP_range) max(vP_range)+0.26]);
ylim([-15 5]);

subplot(2,4,7)

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_SXL(:,1) - vP_range, 1*vP_devs_H2O_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_H2O_SXL(:,2) - vP_range, 1*vP_devs_H2O_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_H2O_SXL(:,3) - vP_range, 1*vP_devs_H2O_SXL(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.26]);
xlabel('True {\itv_p} (x10^{-3})','FontSize',8);
ylim([-5 5]);

subplot(2,4,6)

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_slow(:,1) - vP_range, 1*vP_devs_H2O_slow(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_H2O_slow(:,2) - vP_range, 1*vP_devs_H2O_slow(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_H2O_slow(:,3) - vP_range, 1*vP_devs_H2O_slow(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.26]);
xlabel('True {\itv_p} (x10^{-3})','FontSize',8);
ylim([-5 5]);

subplot(2,4,8)

plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_H2O_slow_SXL(:,1) - vP_range, 1*vP_devs_H2O_slow_SXL(:,1),'LineWidth',1.1,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_H2O_slow_SXL(:,2) - vP_range, 1*vP_devs_H2O_slow_SXL(:,2),'LineWidth',1.1,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_H2O_slow_SXL(:,3) - vP_range, 1*vP_devs_H2O_slow_SXL(:,3),'LineWidth',1.1,'Color',Colour3);
xlim([min(vP_range) max(vP_range)+0.26]);
xlabel('True {\itv_p} (x10^{-3})','FontSize',8);
ylim([-5 5]);

set(gcf, 'units', 'centimeters','Position', [5 5 17.56 10.54]);

annotation(figure(1),'textbox',[0.072 0.948 0.05 0.045],'String','(A)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.279 0.948 0.06 0.045],'String','(B)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.485 0.948 0.06 0.045],'String','(C)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.691 0.948 0.06 0.045],'String','(D)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.072 0.476 0.06 0.045],'String','(E)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.279 0.476 0.06 0.045],'String','(F)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.485 0.476 0.06 0.045],'String','(G)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);
annotation(figure(1),'textbox',[0.691 0.476 0.06 0.045],'String','(H)','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',9);

set(gcf, 'units', 'centimeters','PaperPosition', [0 0 17.56 10.54]);    % can be bigger than screen
print(gcf, 'Figure_2.png', '-dpng','-r1200');
print(gcf, 'Figure_2.eps', '-depsc', '-r1200' );   %save file as eps w/ 1200dpi