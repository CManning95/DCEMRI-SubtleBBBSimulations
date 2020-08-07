clear; close all;

load('PS_devs_Fp.mat')
load('PS_means_Fp.mat')
load('PS_devs_jitter.mat')
load('PS_means_jitter.mat')
load('PS_devs_VFA.mat')
load('PS_means_VFA.mat')
load('PS_range.mat')
load('PS_H2OPatlak.mat')
load('PS_H2OSXL.mat')

Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];
Colour4 = [0.494 0.184 0.556];

figure(1)

subplot(4,3,1)

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

subplot(4,3,2)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_exclude(:,1) - PS_range, 1*PS_devs_Fp_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_Fp_exclude(:,2) - PS_range, 1*PS_devs_Fp_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_exclude(:,3) - PS_range, 1*PS_devs_Fp_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
title('Bolus injection (with exclusion)');
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(4,3,3)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_slow(:,1) - PS_range, 1*PS_devs_Fp_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_Fp_slow(:,2) - PS_range, 1*PS_devs_Fp_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_slow(:,3) - PS_range, 1*PS_devs_Fp_slow(:,3),'LineWidth',1.3,'Color',Colour3);
title('Slow injection');
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(4,3,4)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_fast(:,1) - PS_range, 1*PS_devs_jitter_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_fast(:,2) - PS_range, 1*PS_devs_jitter_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_fast(:,3) - PS_range, 1*PS_devs_jitter_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-2 2]);
legend({'No delay','+ 6 s','+ 12 s'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 9;

subplot(4,3,5)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_exclude(:,1) - PS_range, 1*PS_devs_jitter_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_exclude(:,2) - PS_range, 1*PS_devs_jitter_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_exclude(:,3) - PS_range, 1*PS_devs_jitter_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(4,3,6)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_slow(:,1) - PS_range, 1*PS_devs_jitter_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_slow(:,2) - PS_range, 1*PS_devs_jitter_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_slow(:,3) - PS_range, 1*PS_devs_jitter_slow(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(4,3,7)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_fast(:,1) - PS_range, 1*PS_devs_VFA_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_fast(:,2) - PS_range, 1*PS_devs_VFA_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_fast(:,3) - PS_range, 1*PS_devs_VFA_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-2 2]);
legend({'B1 corrected','T_1_0 and DCE uncorrected (k=1.3)','T_1_0 and DCE uncorrected (k=0.7)'},'Location','best')
legend('boxoff')
xlabel(['True PS (x10^{-4} min^{-1} )']);

ax = gca;
ax.FontSize = 9;

subplot(4,3,8)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_exclude(:,1) - PS_range, 1*PS_devs_VFA_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_exclude(:,2) - PS_range, 1*PS_devs_VFA_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_exclude(:,3) - PS_range, 1*PS_devs_VFA_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
ylim([-2 2]);
xlabel(['True PS (x10^{-4} min^{-1} )']);

ax = gca;
ax.FontSize = 9;

subplot(4,3,9)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_slow(:,1) - PS_range, 1*PS_devs_VFA_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_slow(:,2) - PS_range, 1*PS_devs_VFA_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_slow(:,3) - PS_range, 1*PS_devs_VFA_slow(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
ylim([-2 2]);
xlabel(['True PS (x10^{-4} min^{-1} )']);

ax = gca;
ax.FontSize = 9;

subplot(4,3,10)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2OPatlak_fast(:,1) - PS_range, 1*PS_devs_H2OPatlak_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_H2OPatlak_fast(:,2) - PS_range, 1*PS_devs_H2OPatlak_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_H2OPatlak_fast(:,3) - PS_range, 1*PS_devs_H2OPatlak_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
title('Bolus injection');
xlim([0 max(PS_range)]);
ylim([-2 2]);
legend({'k_{be} = 1.375 s^{-1}','k_{be} = 2.75 s^{-1}','k_{be} = 5.5 s^{-1}'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 9;

subplot(4,3,11)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2OPatlak_exclude(:,1) - PS_range, 1*PS_devs_H2OPatlak_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_H2OPatlak_exclude(:,2) - PS_range, 1*PS_devs_H2OPatlak_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_H2OPatlak_exclude(:,3) - PS_range, 1*PS_devs_H2OPatlak_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
title('Bolus injection (with exclusion)');
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

subplot(4,3,12)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2OPatlak_slow(:,1) - PS_range, 1*PS_devs_H2OPatlak_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_H2OPatlak_slow(:,2) - PS_range, 1*PS_devs_H2OPatlak_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_H2OPatlak_slow(:,3) - PS_range, 1*PS_devs_H2OPatlak_slow(:,3),'LineWidth',1.3,'Color',Colour3);
title('Slow injection');
xlim([0 max(PS_range)]);
ylim([-2 2]);

ax = gca;
ax.FontSize = 9;

% subplot(4,3,10)
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_SXL_fast(:,1) - PS_range, 1*PS_devs_SXL_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_SXL_fast(:,2) - PS_range, 1*PS_devs_SXL_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_SXL_fast(:,3) - PS_range, 1*PS_devs_SXL_fast(:,3),'LineWidth',1.3,'Color',Colour3);
% ylabel('fitted PS error (x10^{-4} min^{-1} )');
% xlabel(['True PS (x10^{-4} min^{-1} )']);
% xlim([0 max(PS_range)]);
% ylim([-2 2]);
%  legend({'k_{be} = 1.375 s^{-1}','k_{be} = 2.75 s^{-1}','k_{be} = 5.5 s^{-1}'},'Location','best')
%  legend('boxoff')
% 
% ax = gca;
% ax.FontSize = 9;
% 
% subplot(4,3,11)
% 
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_SXL_exclude(:,1) - PS_range, 1*PS_devs_SXL_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_SXL_exclude(:,2) - PS_range, 1*PS_devs_SXL_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_SXL_exclude(:,3) - PS_range, 1*PS_devs_SXL_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
% xlim([0 max(PS_range)]);
% xlabel(['True PS (x10^{-4} min^{-1} )']);
% ylim([-2 2]);
% 
% ax = gca;
% ax.FontSize = 9;
% 
% subplot(4,3,12)
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_SXL_slow(:,1) - PS_range, 1*PS_devs_SXL_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_SXL_slow(:,2) - PS_range, 1*PS_devs_SXL_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_SXL_slow(:,3) - PS_range, 1*PS_devs_SXL_slow(:,3),'LineWidth',1.3,'Color',Colour3);
% xlim([0 max(PS_range)]);
% xlabel(['True PS (x10^{-4} min^{-1} )']);
% ylim([-2 2]);
% 
% ax = gca;
% ax.FontSize = 9;

annotation(figure(1),'textbox',[0.064 0.918 0.05 0.045],'String','a. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.357 0.918 0.06 0.045],'String','a. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.636 0.918 0.06 0.045],'String','a. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.064 0.697 0.06 0.045],'String','b. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.357 0.697 0.06 0.045],'String','b. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.636 0.697 0.06 0.045],'String','b. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.064 0.476 0.06 0.045],'String','c. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.357 0.476 0.06 0.045],'String','c. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.636 0.476 0.06 0.045],'String','c. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.064 0.257 0.06 0.045],'String','d. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.357 0.257 0.06 0.045],'String','d. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.636 0.257 0.06 0.045],'String','d. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);

%set(gcf, 'units', 'centimeters','Position', [5 5 25 28]);
set(gcf, 'units', 'centimeters','PaperPosition', [0 0 25 28]);    % can be bigger than screen 
print(gcf, 'PS_fig.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi

% figure(2) % Water exchange plot
% 
% subplot(2,3,1)
% 
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_H2OPatlak_fast(:,1) - PS_range, 1*PS_devs_H2OPatlak_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_H2OPatlak_fast(:,2) - PS_range, 1*PS_devs_H2OPatlak_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_H2OPatlak_fast(:,3) - PS_range, 1*PS_devs_H2OPatlak_fast(:,3),'LineWidth',1.3,'Color',Colour3);
% ylabel('fitted PS error (x10^{-4} min^{-1} )');
% title('Bolus injection');
% xlim([0 max(PS_range)]);
% ylim([-5 5]);
% legend({'k_{be} = 2.5 s^{-1}','k_{be} = 5 s^{-1}','k_{be} = 10 s^{-1}'},'Location','best')
% legend('boxoff')
% 
% ax = gca;
% ax.FontSize = 9;
% 
% subplot(2,3,2)
% 
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_H2OPatlak_exclude(:,1) - PS_range, 1*PS_devs_H2OPatlak_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_H2OPatlak_exclude(:,2) - PS_range, 1*PS_devs_H2OPatlak_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_H2OPatlak_exclude(:,3) - PS_range, 1*PS_devs_H2OPatlak_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
% title('Bolus injection (with exclusion)');
% xlim([0 max(PS_range)]);
% ylim([-2 2]);
% 
% ax = gca;
% ax.FontSize = 9;
% 
% subplot(2,3,3)
% 
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_H2OPatlak_slow(:,1) - PS_range, 1*PS_devs_H2OPatlak_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_H2OPatlak_slow(:,2) - PS_range, 1*PS_devs_H2OPatlak_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_H2OPatlak_slow(:,3) - PS_range, 1*PS_devs_H2OPatlak_slow(:,3),'LineWidth',1.3,'Color',Colour3);
% title('Slow injection');
% xlim([0 max(PS_range)]);
% ylim([-2 2]);
% 
% ax = gca;
% ax.FontSize = 9;
% 
% subplot(2,3,4)
% 
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_SXL_fast(:,1) - PS_range, 1*PS_devs_SXL_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_SXL_fast(:,2) - PS_range, 1*PS_devs_SXL_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_SXL_fast(:,3) - PS_range, 1*PS_devs_SXL_fast(:,3),'LineWidth',1.3,'Color',Colour3);
% ylabel('fitted PS error (x10^{-4} min^{-1} )');
% xlabel(['True PS (x10^{-4} min^{-1} )']);
% xlim([0 max(PS_range)]);
% ylim([-2 2]);
% % legend({'k_{be} = 2.5 s^{-1}','k_{be} = 5 s^{-1}','k_{be} = 10 s^{-1}'},'Location','best')
% % legend('boxoff')
% 
% ax = gca;
% ax.FontSize = 9;
% 
% subplot(2,3,5)
% 
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_SXL_exclude(:,1) - PS_range, 1*PS_devs_SXL_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_SXL_exclude(:,2) - PS_range, 1*PS_devs_SXL_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_SXL_exclude(:,3) - PS_range, 1*PS_devs_SXL_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
% xlim([0 max(PS_range)]);
% xlabel(['True PS (x10^{-4} min^{-1} )']);
% ylim([-2 2]);
% 
% ax = gca;
% ax.FontSize = 9;
% 
% subplot(2,3,6)
% 
% plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
% errorbar(PS_range, PS_means_SXL_slow(:,1) - PS_range, 1*PS_devs_SXL_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
% errorbar(PS_range + 0.03, PS_means_SXL_slow(:,2) - PS_range, 1*PS_devs_SXL_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
% errorbar(PS_range + 0.06, PS_means_SXL_slow(:,3) - PS_range, 1*PS_devs_SXL_slow(:,3),'LineWidth',1.3,'Color',Colour3);
% xlim([0 max(PS_range)]);
% xlabel(['True PS (x10^{-4} min^{-1} )']);
% ylim([-2 2]);
% 
% ax = gca;
% ax.FontSize = 9;
% 
% annotation(figure(2),'textbox',[0.062 0.92 0.05 0.045],'String','a. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
% annotation(figure(2),'textbox',[0.357 0.92 0.06 0.045],'String','a. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
% annotation(figure(2),'textbox',[0.636 0.92 0.06 0.045],'String','a. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
% annotation(figure(2),'textbox',[0.062 0.447 0.06 0.045],'String','b. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
% annotation(figure(2),'textbox',[0.357 0.447 0.06 0.045],'String','b. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
% annotation(figure(2),'textbox',[0.636 0.457 0.06 0.045],'String','b. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
% 
%  set(gcf, 'units', 'centimeters','Position', [5 5 28 17]);
%  
%  set(gcf, 'units', 'centimeters','PaperPosition', [0 0 28 17]);    % can be bigger than screen 
%  print(gcf, 'PS_SXL_fig.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi
