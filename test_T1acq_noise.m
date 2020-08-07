% Test T1 acquisition error
addpath('DCE_Simulation_Functions');
    
[PhysParam,DCESeqParam,SimParam,T1acqParam] = load_default_params;

PhysParam.S0_tissue = 100;
%PhysParam.T10_tissue_s = 0.917; % baseline T1 tissue (NAWM)
PhysParam.T10_tissue_s = 1.212; % baseline T1 tissue (scGM)
T1acqParam.T1_acq_method = 'VFA';
DCESeqParam.FA_error = 1;

T1acqParam.isFit = [1 1 1]; % which acquisitions to fit
T1acqParam.TR_s = [0.0054 0.0054 0.0054]; % repetition times for T1 acqusition
T1acqParam.FA_nom_rads = [2 5 12] *2*(pi/360); % nominal flip angles for T1 acquisition
T1acqParam.FA_true_rads = DCESeqParam.FA_error * T1acqParam.FA_nom_rads; % derive actual flip angles for T1 acquisition
T1acqParam.isIR = [0 0 0]; % indicates which are IR-SPGR
T1acqParam.TI_s = [NaN NaN NaN]; % Inversion times for T1 acquisition (for HIFI)
T1acqParam.PECentre = [NaN NaN NaN]; % indicates time of centre of k-space
T1acqParam.NReadout = [160 160 160]; % number of readout pulses (Siemens - number of slices)
T1acqParam.T1_SNR = 170; % SNR to achieve Lee, Callaghan et al. 2019 0.7% T1 error test-retest in NAWM (1.3% for scGM)
for n = 1:10000
[PhysParam.T1_tissue_meas_s(n,1),temp,T1acqParam.FA_error_meas,temp2] = MeasureT1(PhysParam.S0_tissue,PhysParam.T10_tissue_s,T1acqParam,T1acqParam.T1_acq_method);
end

mean_T1_error = mean(abs(100 - (PhysParam.T10_tissue_s./PhysParam.T1_tissue_meas_s)*100));
disp(['T1 error = ' num2str(mean_T1_error) ' %'])
