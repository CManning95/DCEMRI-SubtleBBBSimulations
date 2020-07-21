% Generates the plasma flow PS simulation data as shown in Manning et al.
% (2020) Slow injection paper
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;
addpath('DCE_Simulation_Functions');
    
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
  SimParam.water_exch_model = '2S1XA'; % water exchange model to generate signals
 % SimParam.SXLfit = 1; % fit enhancements according to SXL method

%% Generate variable flow/injection delay sims
% ranges of PS and vP to test
PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;

N_PS = size(PS_range,1); %range sizes to test
Fp_ranges = [11 8.25 5.5]; % Plasma flow ranges

%% Generate variable Fp PS graphs
 SimParam.InjectionRate = 'fast';
 
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
 SimParam.t_start_s = 0;
 SimParam.NIgnore = max(SimParam.baselineScans);
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
 
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

%% Generate graphs

% Change scale of PS values (for graphing)
PS_range = PS_range * 1e4;
PS_means_Fp_fast = PS_means_Fp_fast * 1e4;
PS_means_Fp_exclude = PS_means_Fp_exclude * 1e4;
PS_means_Fp_slow = PS_means_Fp_slow * 1e4;
PS_devs_Fp_fast = PS_devs_Fp_fast * 1e4;
PS_devs_Fp_exclude = PS_devs_Fp_exclude * 1e4;
PS_devs_Fp_slow = PS_devs_Fp_slow * 1e4;

save('PS_means_Fp','PS_means_Fp_fast','PS_means_Fp_exclude','PS_means_Fp_slow')
save('PS_devs_Fp','PS_devs_Fp_fast','PS_devs_Fp_exclude','PS_devs_Fp_slow')

Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];

subplot(1,3,1)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_fast(:,1) - PS_range, 1*PS_devs_Fp_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_Fp_fast(:,2) - PS_range, 1*PS_devs_Fp_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_fast(:,3) - PS_range, 1*PS_devs_Fp_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
title('Bolus injection');
xlim([0 max(PS_range)]);
ylim([-2 2]);
legend({'F_p = 11 ml 100g^{-1}min^{-1}','F_p = 8.25 ml 100g^{-1}min^{-1}','F_p = 5.5 ml 100g^{-1}min^{-1}'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 9;

subplot(1,3,2)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_exclude(:,1) - PS_range, 1*PS_devs_Fp_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_Fp_exclude(:,2) - PS_range, 1*PS_devs_Fp_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_exclude(:,3) - PS_range, 1*PS_devs_Fp_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
title('Bolus injection (with exclusion)');
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(1,3,3)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_slow(:,1) - PS_range, 1*PS_devs_Fp_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_Fp_slow(:,2) - PS_range, 1*PS_devs_Fp_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_slow(:,3) - PS_range, 1*PS_devs_Fp_slow(:,3),'LineWidth',1.3,'Color',Colour3);
title('Slow injection');
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

set(gcf, 'units', 'centimeters','Position', [5 5 28 10]);
