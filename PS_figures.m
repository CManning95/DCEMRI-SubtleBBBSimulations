clear; close all;

load('PS_devs.mat')
load('PS_means.mat')
load('PS_range.mat')

Colour1  = [0 0.447 0.741 0.5];
Colour2 = [0.85 0.325 0.098 0.5];
Colour3 = [0.929 0.694 0.125 0.5];
Colour4 = [0.494 0.184 0.556];

figure(1)

subplot(5,3,1)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_fast(:,1) - PS_range, 1*PS_devs_Fp_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_Fp_fast(:,2) - PS_range, 1*PS_devs_Fp_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_fast(:,3) - PS_range, 1*PS_devs_Fp_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
title('Bolus injection');
xlim([0 max(PS_range)]);
ylim([-3.2 3.2]);
legend({'F_p = 18 ml 100g^{-1}min^{-1}','F_p = 9 ml 100g^{-1}min^{-1}','F_p = 4.5 ml 100g^{-1}min^{-1}'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 10;

subplot(5,3,2)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_exclude(:,1) - PS_range, 1*PS_devs_Fp_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_Fp_exclude(:,2) - PS_range, 1*PS_devs_Fp_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_exclude(:,3) - PS_range, 1*PS_devs_Fp_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
title('Bolus injection (with exclusion)');
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,3)

plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_Fp_slow(:,1) - PS_range, 1*PS_devs_Fp_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_Fp_slow(:,2) - PS_range, 1*PS_devs_Fp_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_Fp_slow(:,3) - PS_range, 1*PS_devs_Fp_slow(:,3),'LineWidth',1.3,'Color',Colour3);
title('Slow injection');
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,4)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_fast(:,1) - PS_range, 1*PS_devs_jitter_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_fast(:,2) - PS_range, 1*PS_devs_jitter_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_fast(:,3) - PS_range, 1*PS_devs_jitter_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-3.2 3.2]);
legend({'No delay','+ 6 s','+ 12 s'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 10;

subplot(5,3,5)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_exclude(:,1) - PS_range, 1*PS_devs_jitter_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_exclude(:,2) - PS_range, 1*PS_devs_jitter_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_exclude(:,3) - PS_range, 1*PS_devs_jitter_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,6)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_jitter_slow(:,1) - PS_range, 1*PS_devs_jitter_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_jitter_slow(:,2) - PS_range, 1*PS_devs_jitter_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_jitter_slow(:,3) - PS_range, 1*PS_devs_jitter_slow(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,7)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_fast(:,1) - PS_range, 1*PS_devs_H2O_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_H2O_fast(:,2) - PS_range, 1*PS_devs_H2O_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_H2O_fast(:,3) - PS_range, 1*PS_devs_H2O_fast(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
errorbar(PS_range + 0.09, PS_means_H2O_fast(:,4) - PS_range, 1*PS_devs_H2O_fast(:,4),'LineWidth',1.3,'Color',Colour4);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-3.2 3.2]);
legend({'FXL','k_{be} = 2.5 s^{-1}','k_{be} = 5 s^{-1}','k_{be} = 10 s^{-1}'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 10;

subplot(5,3,8)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_exclude(:,1) - PS_range, 1*PS_devs_H2O_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_H2O_exclude(:,2) - PS_range, 1*PS_devs_H2O_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_H2O_exclude(:,3) - PS_range, 1*PS_devs_H2O_exclude(:,3),'LineWidth',1.3,'Color',Colour3); hold on;
errorbar(PS_range + 0.09, PS_means_H2O_exclude(:,4) - PS_range, 1*PS_devs_H2O_exclude(:,4),'LineWidth',1.3,'Color',Colour4);
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,9)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_H2O_slow(:,1) - PS_range, 1*PS_devs_H2O_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_H2O_slow(:,2) - PS_range, 1*PS_devs_H2O_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_H2O_slow(:,3) - PS_range, 1*PS_devs_H2O_slow(:,3),'LineWidth',1.3,'Color',Colour3);
errorbar(PS_range + 0.09, PS_means_H2O_slow(:,4) - PS_range, 1*PS_devs_H2O_slow(:,4),'LineWidth',1.3,'Color',Colour4);
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,10)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_fast(:,1) - PS_range, 1*PS_devs_VFA_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_fast(:,2) - PS_range, 1*PS_devs_VFA_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_fast(:,3) - PS_range, 1*PS_devs_VFA_fast(:,3),'LineWidth',1.3,'Color',Colour3);
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-3.2 3.2]);
legend({'T_1_0 and DCE corrected','T_1_0 corrected only','T_1_0 and DCE uncorrected'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 10;

subplot(5,3,11)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_exclude(:,1) - PS_range, 1*PS_devs_VFA_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_exclude(:,2) - PS_range, 1*PS_devs_VFA_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_exclude(:,3) - PS_range, 1*PS_devs_VFA_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,12)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_VFA_slow(:,1) - PS_range, 1*PS_devs_VFA_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_VFA_slow(:,2) - PS_range, 1*PS_devs_VFA_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_VFA_slow(:,3) - PS_range, 1*PS_devs_VFA_slow(:,3),'LineWidth',1.3,'Color',Colour3);
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,13)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_T1assumed_fast(:,1) - PS_range, 1*PS_devs_T1assumed_fast(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_T1assumed_fast(:,2) - PS_range, 1*PS_devs_T1assumed_fast(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_T1assumed_fast(:,3) - PS_range, 1*PS_devs_T1assumed_fast(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True PS (x10^{-4} min^{-1} )');
ylabel('fitted PS error (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-3.2 3.2]);
legend({'T_{10} accurate','T_{10} error = -20%','T_{10} error = +20%'},'Location','best')
legend('boxoff')

ax = gca;
ax.FontSize = 10;

subplot(5,3,14)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_T1assumed_exclude(:,1) - PS_range, 1*PS_devs_T1assumed_exclude(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_T1assumed_exclude(:,2) - PS_range, 1*PS_devs_T1assumed_exclude(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_T1assumed_exclude(:,3) - PS_range, 1*PS_devs_T1assumed_exclude(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True PS (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

subplot(5,3,15)
plot(PS_range,zeros(size(PS_range)),'k:','DisplayName','True PS','HandleVisibility','off'); hold on;
errorbar(PS_range, PS_means_T1assumed_slow(:,1) - PS_range, 1*PS_devs_T1assumed_slow(:,1),'LineWidth',1.3,'Color',Colour1); hold on;
errorbar(PS_range + 0.03, PS_means_T1assumed_slow(:,2) - PS_range, 1*PS_devs_T1assumed_slow(:,2),'LineWidth',1.3,'Color',Colour2); hold on;
errorbar(PS_range + 0.06, PS_means_T1assumed_slow(:,3) - PS_range, 1*PS_devs_T1assumed_slow(:,3),'LineWidth',1.3,'Color',Colour3);
xlabel('True PS (x10^{-4} min^{-1} )');
xlim([0 max(PS_range)]);
ylim([-1.6 1.6]);

ax = gca;
ax.FontSize = 10;

annotation(figure(1),'textbox',[0.07 0.9 0.05 0.045],'String','a. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.362 0.9 0.06 0.045],'String','a. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.643 0.9 0.06 0.045],'String','a. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.07 0.729 0.06 0.045],'String','b. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.362 0.729 0.06 0.045],'String','b. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.643 0.729 0.06 0.045],'String','b. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.07 0.556 0.06 0.045],'String','c. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.362 0.556 0.06 0.045],'String','c. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.643 0.556 0.06 0.045],'String','c. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.07 0.383 0.06 0.045],'String','d. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.362 0.383 0.06 0.045],'String','d. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.643 0.383 0.06 0.045],'String','d. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.07 0.21 0.06 0.045],'String','e. i','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.362 0.21 0.06 0.045],'String','e. ii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);
annotation(figure(1),'textbox',[0.643 0.21 0.06 0.045],'String','e. iii','LineStyle','none','FitBoxToText','off','fontweight','bold','FontSize',13);



set(gcf, 'PaperPosition', [0 0 28 40]);    % can be bigger than screen 
print(gcf, 'PS_fig.png', '-dpng', '-r300' );   %save file as PNG w/ 300dpi
