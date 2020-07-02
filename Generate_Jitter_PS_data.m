% Generates the injection delay PS simulation data as shown in Manning et al.
% (2020) Slow injection paper
% This is purely for aesthetics - all simulations can be run from the GUI
% for ease of use

clear; close all;

addpath('DCE_Simulation_Functions');

[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;
SimParam.water_exch_model = '2S1XA'; % water exchange model to generate signals
SimParam.SXLfit = 0; % fit enhancements according to SXL method

PS_range = linspace(SimParam.min_PS,SimParam.max_PS,10)'+1e-8;
vP_fixed = PhysParam.vP_fixed;
N_PS = size(PS_range,1); %range sizes to test
Delay_ranges = [0 10 20]; % Injection delay ranges

%% Generate variable jitter PS - fast injection
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

%% Generate variable jitter PS - fast injection (exclude)
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
 
 %% Generate variable jitter PS - slow injection
 SimParam.t_start_s = 0;
 SimParam.InjectionRate = 'slow';
 SimParam.baselineScans = [1:3]; % datapoints to use for calculating base signal
 SimParam.NIgnore = max(SimParam.baselineScans);
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
 
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

%% convert units, save results
PS_range = PS_range * 1e4;
PS_means_jitter_fast = PS_means_jitter_fast * 1e4;
PS_means_jitter_exclude = PS_means_jitter_exclude * 1e4;
PS_means_jitter_slow = PS_means_jitter_slow * 1e4;
PS_devs_jitter_fast = PS_devs_jitter_fast * 1e4;
PS_devs_jitter_exclude = PS_devs_jitter_exclude * 1e4;
PS_devs_jitter_slow = PS_devs_jitter_slow * 1e4;

save('PS_means_jitter','PS_means_jitter_fast','PS_means_jitter_exclude','PS_means_jitter_slow')
save('PS_devs_jitter','PS_devs_jitter_fast','PS_devs_jitter_exclude','PS_devs_jitter_slow')

%% Plot figure of results
Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];

figure()
subplot(1,3,1)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_fast(:,1) - PS_range, 1*PS_devs_jitter_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_fast(:,2) - PS_range, 1*PS_devs_jitter_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_fast(:,3) - PS_range, 1*PS_devs_jitter_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
title('Bolus injection');
xlim([0 max(PS_range)]);
ylim([-4 4]);
legend({'No delay','+ 2 s','+ 4 s'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 9;

subplot(1,3,2)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_exclude(:,1) - PS_range, 1*PS_devs_jitter_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_exclude(:,2) - PS_range, 1*PS_devs_jitter_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_exclude(:,3) - PS_range, 1*PS_devs_jitter_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
title('Bolus injection (with exclusion)');
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(1,3,3)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_slow(:,1) - PS_range, 1*PS_devs_jitter_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_slow(:,2) - PS_range, 1*PS_devs_jitter_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_slow(:,3) - PS_range, 1*PS_devs_jitter_slow(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
title('Slow injection');
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

set(gcf, 'units', 'centimeters','Position', [5 5 28 10]);