function[vP_fit,PS_fit_perMin]=test_master_single_sim(PhysParam,SeqParam,SimParam)
% Master single simulation for a ROI to assess Patlak for BBB permeability
% vP_fit = fitted vP
% PS_fit_perMin = fitted PS

%% Derive additional parameters
vB = PhysParam.vP/(1-PhysParam.Hct); % volume fraction for blood
vI = 1 - (vB + PhysParam.vE); % volume fraction that is intracellular space (assuming vB + vE + vI = 1)
V = [vB PhysParam.vE vI].'; % volume fractions for blood, EES, cell
p = V; % array of water population fractions (blood, EES, cell), assuming equal density (i.e. partition coefficient not taken into consideration)

% calculate initial equilibrium magnetisation vector, i.e. [M0(blood) M0(EES) M0(intracellular)]
M0 = PhysParam.S0_tissue * p; % assume magnetisation is proportional to the population fraction, scaled to give the specified total equilibrium signal S0

% For water exchange models, calculate compartmental T1 values by assuming
% FXL is valid pre-contrast so that R10 is a weighted average of
% compartment values. We also assume that intra- and extra-cellular R10 are
% the same
T10_ES_s = (M0(2) + M0(3)) / ( ((1/PhysParam.T10_tissue_s)*(M0(1)+M0(2)+M0(3))) - ( (1/PhysParam.T10_blood_s)*M0(1)) ); %calculate T10 for combined extravascular space
T10_s = [PhysParam.T10_blood_s T10_ES_s T10_ES_s]; % vector of compartmental T1s for blood, EES, intracellular; assume T10 is equal in EES and cells (equivalent to assuming FX between interstitium and cells)

% If SXLfit not selected, default to traditional Patlak analysis
if isfield(SimParam,'SXLfit') == 0; 
    SimParam.SXLfit = 0;
end

if isfield(SimParam,'Plot_extra_figs') == 0; 
    SimParam.Plot_extra_figs = 0;
end

%% Simulating synthetic AIF and tissue signals
% start by generating AIF concentration curve
switch SimParam.InjectionRate
    case 'fast' % Use Parker model to generate synthetic 'fast' injection
        [timepoints_full_s, Cp_AIF_mM] = DCEFunc_getParkerModAIF(SimParam.t_res_full_s, SeqParam.t_acq_s, SimParam.t_start_s, PhysParam.Hct, 'OG Parker-MSS3'); % Generate AIF from modified Parker model
    case 'slow' % use slow injection AIF from patient data to generate 'slow' injection
        t_AIF_input_s = ((1:1:SimParam.InputAIFDCENFrames) * (SimParam.tRes_InputAIF_s)) - (SimParam.tRes_InputAIF_s/2).'; % Assume concentration is measured at centre of acquisition
        initial_timepoints_full_s = (0:SimParam.t_res_full_s:(SeqParam.t_acq_s - SimParam.t_start_s)).'; % full range of upsampled time points (subtract additional pre-injection delay from end of AIF)
        Cp_AIF_mM = interp1(t_AIF_input_s,SimParam.Cp_AIF_mM,initial_timepoints_full_s,'spline','extrap'); % interpolate AIF to high temporal resolution
        Cp_AIF_mM(1:(3/SimParam.InputAIFDCENFrames)*size(initial_timepoints_full_s)) = 0; %Set all pre-contrast points to zero (as interpolation produces non-zero pre contrast values)
        % adding zeroes to start for injection delay
        t_preContrast_s=fliplr(SimParam.t_start_s-SimParam.t_res_full_s/2:-SimParam.t_res_full_s:0).'; %pre-injection time points (calculate backwards to zero, then reverse, so that time interval is constant)
        Cp_AIF_mM=[zeros(size(t_preContrast_s)) ; Cp_AIF_mM]; % add pre-contrast concentration values (zeros)
        timepoints_full_s=[ t_preContrast_s ; initial_timepoints_full_s + SimParam.t_start_s]; %add pre-contrast time points to time values
        
end

% Generate 'measured' AIF for DCE fitting (add temporal delay to AIF used
% to generate tissue concentrations)
AIF_venous_delay_s = fliplr(SimParam.venous_delay_s-SimParam.t_res_full_s/2:-SimParam.t_res_full_s:0).'; % flip and add zeros at start for venous delay
meas_AIF_mM = [zeros(size(AIF_venous_delay_s)) ; Cp_AIF_mM]; % add delayed concentration values (zeros)
meas_AIF_mM = meas_AIF_mM(1:end - size(AIF_venous_delay_s),:); % cuts off end to make same size as Cp_AIF_mM 

%% generate tissue concentration curves with 2CXM; c_cp_mM = capillary plasma local conc, c_e_mM = EES local conc
[Ct_mM, IRF, c_cp_mM, c_ees_mM] = DCEFunc_PKP2Conc(SimParam.t_res_full_s, Cp_AIF_mM, struct('vP',PhysParam.vP,'vE',PhysParam.vE,'PS_perMin', PhysParam.PS_perMin, 'FP_mlPer100gPerMin', PhysParam.FP_mlPer100gPerMin),SimParam.syn_model,[]);
Ct_cp_mM=PhysParam.vP*c_cp_mM; %capillary plasma contribution to overall Ct
Ct_ees_mM=PhysParam.vE*c_ees_mM; %EES contribution to overall Ct
c_ev_mM=PhysParam.vE*c_ees_mM / (PhysParam.vE+vI); %total extravascular Gd local concentration (i.e. quantity normalised to EES + intracellular volume)

%% Generate synthetic AIF/tissue signals
drift = ((SimParam.drift_pctPerMin/100)/(60))*timepoints_full_s; %fractional signal change due to drift at each synthetic data time point, relative to time zero

% generate synthetic 'measured' AIF signal at full set of timepoints (assuming fast water exchange within blood - i.e. between red blood cells/plasma)
enh_AIF_pct = DCEFunc_Conc2Enh_SPGR((1-PhysParam.Hct)*(meas_AIF_mM),PhysParam.T10_blood_s, SeqParam.TR_s, SeqParam.TE_s, SeqParam.blood_FA_true_deg, SeqParam.r1_per_mM_per_s, SeqParam.r2_per_mM_per_s); % percentage enhancement of AIF signal
SI_AIF = DCEFunc_getSPGRSignal(PhysParam.S0_blood, PhysParam.T10_blood_s, PhysParam.T2s0_blood_s, SeqParam.TR_s, SeqParam.TE_s, SeqParam.blood_FA_true_deg) * (1+enh_AIF_pct/100); % AIF signal
SI_AIF_drifted = SI_AIF.*(1+drift); % apply drift to AIF signal

% generate synthetic tissue signal at full set of timepoints 
% NB we need to use the local capillary concentration c_cp_mM here, NOT the AIF concentration Cp_AIF_mM
enh_tissue_pct = H2OEx_Conc2Enh_SPGR(... % first calculate enhancement based on chosen compartmental relaxation model
    (1-PhysParam.Hct)*c_cp_mM,c_ees_mM,V,p,T10_s.',PhysParam.T2s0_tissue_s,SeqParam.TR_s,SeqParam.TE_s,SeqParam.tissue_FA_true_deg,PhysParam.kbe_perS,PhysParam.kie_perS,SeqParam.r1_per_mM_per_s,SeqParam.r2_per_mM_per_s,SimParam.water_exch_model); % percentage enhancement in tissue signal
SI_tissue = H2OEx_getSPGRSignal(M0,T10_s.',PhysParam.T2s0_tissue_s,SeqParam.TR_s,SeqParam.TE_s,SeqParam.tissue_FA_true_deg,PhysParam.kbe_perS,PhysParam.kie_perS,SimParam.water_exch_model) * (1+enh_tissue_pct/100); % tissue signal
SI_tissue_drifted = SI_tissue.*(1+drift); % apply drift to tissue signal

%% Downsample synthetic signal to generate measured signal
[t_sample, SI_AIF_sample] = DCEFunc_downSample(timepoints_full_s,SI_AIF_drifted,SeqParam.t_res_sample_s,'nearestPoint');
[temp,SI_tissue_sample] = DCEFunc_downSample(timepoints_full_s,SI_tissue_drifted,SeqParam.t_res_sample_s,'nearestPoint');

%% add noise to downsampled tissue signal
SI_AIF_sample = repmat(SI_AIF_sample,1,SimParam.N_repetitions);
SI_tissue_sample = repmat(SI_tissue_sample,1,SimParam.N_repetitions); %generate multiple copies so that noise effects can be simulated
sigma_signal_noise = SI_tissue(1)/SimParam.SNR; % standard deviation of noise
tissue_signal_noise = sigma_signal_noise * randn(size(SI_tissue_sample)); % random array of noise, same size as SI_tissue_sample
SI_tissue_sample_noisy = SI_tissue_sample + tissue_signal_noise; % add noise to sampled signal

%% Calculate sampled AIF and tissue enhancements from sampled signal
enh_AIF_sample_pct = DCEFunc_Sig2Enh(SI_AIF_sample,SimParam.baselineScans); % sampled AIF signal enhancement percentages
enh_tissue_sample_pct = DCEFunc_Sig2Enh(SI_tissue_sample_noisy,SimParam.baselineScans); % sampled tissue signal enhancement percentages

%% Sort measured T1 and FA input for concentration calculations
if size(PhysParam.T1_tissue_meas_s) == 1; % if using Accurate/Assumed T1 acq, use same T1 for each iteration
    T1_blood_Enh2Conc_input = ones(1,SimParam.N_repetitions)*PhysParam.T1_blood_meas_s;
    T1_tissue_Enh2Conc_input = ones(1,SimParam.N_repetitions)*PhysParam.T1_tissue_meas_s;   
else % if using VFA T1 acq, use different T1 for each iteration (it will be of size N_repetitions anyway)
    T1_blood_Enh2Conc_input = PhysParam.T1_blood_meas_s';
    T1_tissue_Enh2Conc_input = PhysParam.T1_tissue_meas_s';
end
if size(SeqParam.blood_FA_meas_deg) == 1;
    FA_blood_meas_input = (ones(1,SimParam.N_repetitions)*SeqParam.blood_FA_meas_deg);
    FA_tissue_meas_input = (ones(1,SimParam.N_repetitions)*SeqParam.tissue_FA_meas_deg);
else
    FA_blood_meas_input = SeqParam.blood_FA_meas_deg';
    FA_tissue_meas_input = SeqParam.tissue_FA_meas_deg';
end
    
%% Calculate sampled AIF and tissue concentrations using enhancements at sampling points (note we use "measured" FA and T1 values to allow for possible errors)
if SeqParam.r2_per_mM_per_s == 0;
    mode = 'analytical';
else
    mode = 'numeric';
end
Cp_AIF_sample_mM = ((1/(1-PhysParam.Hct))*DCEFunc_Enh2Conc_SPGR(enh_AIF_sample_pct,T1_blood_Enh2Conc_input,SeqParam.TR_s,SeqParam.TE_s,FA_blood_meas_input,SeqParam.r1_per_mM_per_s,SeqParam.r2_per_mM_per_s,mode));

%% Fit enhancement curves using SXL fitting method or Patlak model
if SimParam.SXLfit == 1;
    opts.init_vP = 0.0001; % initial vP value to test for SXL fitting
    opts.init_PS_perMin = 1e-8; % initial PS value to test for SXL fitting
    opts.NIgnore = SimParam.NIgnore;
    [PatlakResults, enhModelFit_pct] = DCEFunc_fitPatlak_waterEx(SeqParam.t_res_sample_s,enh_tissue_sample_pct,Cp_AIF_sample_mM,PhysParam.Hct,T1_tissue_Enh2Conc_input,T1_blood_Enh2Conc_input,SeqParam.TR_s,SeqParam.TE_s,FA_tissue_meas_input,SeqParam.r1_per_mM_per_s,SeqParam.r2_per_mM_per_s,opts);
    Ct_sample_mM = PatlakResults.Ct_SXL_mM;
    Ct_fit_mM = mean(Ct_sample_mM,2);
    
else
    Ct_sample_mM = DCEFunc_Enh2Conc_SPGR(enh_tissue_sample_pct,T1_tissue_Enh2Conc_input,SeqParam.TR_s,SeqParam.TE_s,FA_tissue_meas_input,SeqParam.r1_per_mM_per_s,SeqParam.r2_per_mM_per_s,mode); 
    [PatlakResults, Ct_fit_mM] = DCEFunc_fitModel(SeqParam.t_res_sample_s,Ct_sample_mM,Cp_AIF_sample_mM,'PatlakFastMultiAIF',struct('NIgnore',SimParam.NIgnore));
end

%% plot results
if SimParam.Plot_extra_figs == 1;
    figure(1); set(gcf,'Units','normalized','outerposition',[0 0 1 1]);
    suptitle(['v_p=' num2str(PhysParam.vP) ', PS (per Min)=' num2str(PhysParam.PS_perMin) ', Fp = ' num2str(PhysParam.FP_mlPer100gPerMin) ', model = ' SimParam.water_exch_model]);
    subplot(3,2,1)
    plot(timepoints_full_s,meas_AIF_mM,'k-',timepoints_full_s,Cp_AIF_mM,'r-',timepoints_full_s,c_cp_mM,'m:'); hold on;
    errorbar(t_sample,mean(Cp_AIF_sample_mM,2),1*std(Cp_AIF_sample_mM,0,2),'.');
    title(['AIF Conc (\deltat=' num2str(SimParam.t_res_full_s) ')']);
    ylabel('C_p (AIF) / mMol');
    legend({'VIF', 'AIF','sampled VIF','c_c_p'}, 'Location', 'northeastoutside','Fontsize',10)
    subplot(3,2,3)
    plot(timepoints_full_s,enh_AIF_pct,'k-'); hold on;
    errorbar(t_sample,mean(enh_AIF_sample_pct,2),1*std(enh_AIF_sample_pct,0,2),'.');
    title('VIF enhancement');
    ylabel('enh / %');
    legend({'VIF enh %','sampled VIF enh %'},'Location','northeastoutside','Fontsize',10)
    subplot(3,2,2)
    plot(timepoints_full_s,Ct_mM,'k-',timepoints_full_s,Ct_cp_mM,'r--',timepoints_full_s,Ct_ees_mM,'g--'); hold on
    errorbar(t_sample,mean(Ct_sample_mM,2),1*std(Ct_sample_mM,0,2),'.');
    title(['Tissue Conc contributions']);
    ylabel('C_t / mMol');
    legend({'Full C_t','vascular conc', 'extravascular conc', 'Sampled C_t'}, 'Location', 'northeastoutside', 'Fontsize',10)
    subplot(3,2,4)
    title(['Ccp, EES and EV conc']);
    yyaxis left
    plot(timepoints_full_s,c_cp_mM,'r-');
    ylabel('c_c_p / mMol');
    ylim([0 10]);
    yyaxis right
    plot(timepoints_full_s,c_ees_mM,'g-',timepoints_full_s,c_ev_mM,'b-');
    ylabel('c_e_e_s and c_e_v / mMol');
    ylim([0 0.04]);
    legend({'c_c_p','c_e_e_s','c_e_v'}, 'Location', 'northeastoutside','Fontsize',10)
    subplot(3,2,6)
    plot(timepoints_full_s,enh_tissue_pct,'k:'); hold on;
    if SimParam.SXLfit == 1;
        errorbar(t_sample,mean(enhModelFit_pct,2),1*std(enhModelFit_pct,0,2),'.');
    else
        errorbar(t_sample,mean(enh_tissue_sample_pct,2),1*std(enh_tissue_sample_pct,0,2),'.');
    end
    title('tissue enhancement');
    ylabel('enh tissue / %');
    xlabel('time (s)')
    legend({'full tissue enh %','sampled tissue enh %'}, 'Location','northeastoutside','Fontsize',10)
    subplot(3,2,5)
    plot(timepoints_full_s,Ct_mM,'k:'); hold on
    errorbar(t_sample,mean(Ct_sample_mM,2),1*std(Ct_sample_mM,0,2),'.'); hold on;
    plot(t_sample,mean(Ct_fit_mM,2));
    legend({'full C_t','measured C_t','fitted C_t (Patlak)'},'Location','northeastoutside','Fontsize',10);
    title(['Patlak: mean v_p=' num2str(mean(PatlakResults.vP)) ', mean PS(per Min)=' num2str(mean(PatlakResults.PS_perMin))]);
    ylabel('C_t / mMol');
    xlabel('time (s)');
    %ylim([min(Ct_sample_mM(:,1)) max(Ct_sample_mM(:,1))]);
    pause(0.01)

end
vP_fit = PatlakResults.vP;
PS_fit_perMin = PatlakResults.PS_perMin;

end