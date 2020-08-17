% Generates the injection delay PS simulation data as shown in Manning et al.
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

N_PS = size(PS_range,1); %range sizes to test
N_vP = size(vP_range,1); %range sizes to test
Delay_ranges = [0 4 8 12]; % Injection delay ranges

% T1 acquisition
[PhysParam.T1_blood_meas_s,temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_blood,PhysParam.T10_blood_s,T1acqParam,T1acqParam.T1_acq_method);
for m = 1:N_PS
    for n = 1:SimParam.N_repetitions
        [T1_tissue_meas_s(n,m),temp,temp2,temp3] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
    end
end

%% Generate variable jitter PS - fast injection
for i = 1:size(Delay_ranges,2);
     SimParam.t_start_s = 119 + Delay_ranges(i);
     for i_PS = 1:N_PS
         PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
         PhysParam.vP = vP_fixed(1);
         PhysParam.PS_perMin = PS_range(i_PS);
         [temp, PS_fit_jitter_fast(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
     end
     for i_vP = 1:N_vP
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
 
 %% Generate variable jitter PS - slow injection
 SimParam.t_start_s = 0;
 SimParam.InjectionRate = 'slow';
 
 load('Slow_Cp_AIF_mM.mat') % load example slow injection VIF
 SimParam.Cp_AIF_mM = Cp_AIF_mM;
 SimParam.tRes_InputAIF_s = 39.62; % original time resolution of AIFs
 SimParam.InputAIFDCENFrames = 32; % number of time points
 
   for i = 1:size(Delay_ranges,2);
      SimParam.t_start_s = Delay_ranges(i);
      for i_PS = 1:N_PS
          PhysParam.T1_tissue_meas_s = T1_tissue_meas_s(:,i_PS);
          PhysParam.vP = vP_fixed(1);
          PhysParam.PS_perMin = PS_range(i_PS);
          [temp, PS_fit_jitter_slow(:,i_PS)] = master_single_sim(PhysParam,DCESeqParam,SimParam);
      end
      for i_vP = 1:N_vP
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

%% convert units, save results
PS_range = PS_range * 1e4;
PS_means_jitter_fast = PS_means_jitter_fast * 1e4;
PS_means_jitter_slow = PS_means_jitter_slow * 1e4;
PS_devs_jitter_fast = PS_devs_jitter_fast * 1e4;
PS_devs_jitter_slow = PS_devs_jitter_slow * 1e4;

vP_range = vP_range * 1e3;
vP_means_jitter_fast = vP_means_jitter_fast * 1e3;
vP_means_jitter_slow = vP_means_jitter_slow * 1e3;
vP_devs_jitter_fast = vP_devs_jitter_fast * 1e3;
vP_devs_jitter_slow = vP_devs_jitter_slow * 1e3;

save('PS_means_jitter','PS_means_jitter_fast','PS_means_jitter_slow')
save('PS_devs_jitter','PS_devs_jitter_fast','PS_devs_jitter_slow')
save('vP_means_jitter','vP_means_jitter_fast','vP_means_jitter_slow')
save('vP_devs_jitter','vP_devs_jitter_fast','vP_devs_jitter_slow')

%% Plot figure of results
Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];
Colour4 = [0.4940 0.1840 0.5560 0.5];

figure()
subplot(2,2,1)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_fast(:,1) - PS_range, 1*PS_devs_jitter_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_fast(:,2) - PS_range, 1*PS_devs_jitter_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_fast(:,3) - PS_range, 1*PS_devs_jitter_fast(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_fast(:,4) - PS_range, 1*PS_devs_jitter_fast(:,4),'LineWidth',1.3,'Color',Colour4);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlabel(['True PS (x10^{-4} min^{-1} )']);
title('Bolus injection');
xlim([0 max(PS_range)]);
ylim([-2 2]);
legend({'No delay','+ 4 s','+ 8 s','+ 12 s'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 9;

subplot(2,2,2)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_slow(:,1) - PS_range, 1*PS_devs_jitter_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_slow(:,2) - PS_range, 1*PS_devs_jitter_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.12, PS_means_jitter_slow(:,3) - PS_range, 1*PS_devs_jitter_slow(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
errorbar(PS_range + 0.18, PS_means_jitter_slow(:,4) - PS_range, 1*PS_devs_jitter_slow(:,4),'LineWidth',1.3,'Color',Colour4);
xlim([0 max(PS_range)]);
title('Slow injection');
xlabel(['True PS (x10^{-4} min^{-1} )']);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(2,2,3)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_fast(:,1) - vP_range, 1*vP_devs_jitter_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_jitter_fast(:,2) - vP_range, 1*vP_devs_jitter_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_jitter_fast(:,3) - vP_range, 1*vP_devs_jitter_fast(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
errorbar(vP_range + 0.39, vP_means_jitter_fast(:,4) - vP_range, 1*vP_devs_jitter_fast(:,4),'LineWidth',1.3,'Color',Colour4);
ylabel('fitted v_p error (x10^{-3})');
xlabel(['True v_p (x10^{-3})']);
xlim([min(vP_range) max(vP_range)+0.26]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(2,2,4)
plot(vP_range,zeros(size(vP_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(vP_range, vP_means_jitter_slow(:,1) - vP_range, 1*vP_devs_jitter_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(vP_range + 0.13, vP_means_jitter_slow(:,2) - vP_range, 1*vP_devs_jitter_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(vP_range + 0.26, vP_means_jitter_slow(:,3) - vP_range, 1*vP_devs_jitter_slow(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
errorbar(vP_range + 0.39, vP_means_jitter_slow(:,4) - vP_range, 1*vP_devs_jitter_slow(:,4),'LineWidth',1.3,'Color',Colour4);
xlim([min(vP_range) max(vP_range)+0.26]);
xlabel(['True v_p (x10^{-3})']);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

annotation(figure(1),'textbox',[0.064 0.935 0.07 0.045],'String','a. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.488 0.935 0.07 0.045],'String','a. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.064 0.460 0.07 0.045],'String','b. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
annotation(figure(1),'textbox',[0.488 0.460 0.07 0.045],'String','b. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',11);
set(gcf, 'units', 'centimeters','Position', [5 5 18 15]);
set(gcf, 'units', 'centimeters','PaperPosition', [0 0 15 15]);    % can be bigger than screen 
print(gcf, 'Jitter_figure.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi