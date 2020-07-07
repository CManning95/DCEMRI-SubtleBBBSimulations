function[T1_meas,S0_meas,k_meas,modelFit] = MeasureT1(S0,T1,acqParam,T1_acq_method,varargin)
% MeasureT1 simulates T1 acquisition using either Variable Flip Angle (VFA)
% or HIFI methods to measure T1
% Input Parameters:
% S0 - base signal
% T1 - actual base T1 to measure
% --variable input --
% assumed_T1 - assumed T1 in s, only for 'Assumed' acquisition method

% acqParam: struct containing T1 measurement acquisition parameters:
%   TR_s - array of repetition times
%   FA_acq_true_rads - array of flip actual flip angles in radians
%   isIR - array to select which are IR
%   isFit - array to select which to fit
%   FA_acq_nom_rads - array of nominal flip angles in radians
%   TI_s - array of inversion times (if IR)
%   PECentre - array indicating time of centre of k-space 
%   NReadout - array of number of readout pulses
%   NTry - Attempts at fitting

%T1_acq_method - Acquisition method - VFA or HIFI

switch T1_acq_method
    case 'Accurate' % Assume T1 is known exactly
        T1_meas = T1;
        S0_meas = NaN;
        k_meas = NaN;
        modelFit = NaN;
    case {'VFA','HIFI'}
        NIR=sum(acqParam.isFit,2);
        NT1=size(acqParam.isFit,2);
        % generate T1 scan intensities
        SI=nan(1,NT1); %pre-allocate
        T1_noise=nan(1,NT1);
        
        SI(~acqParam.isIR) = SPGRFormula(S0,T1,acqParam.TR_s(~acqParam.isIR),acqParam.FA_true_rads(~acqParam.isIR)); % signal of SPGR components
        SI(acqParam.isIR==1) = deichmannFormula(S0,T1,acqParam.TR_s(acqParam.isIR==1),acqParam.TI_s(acqParam.isIR==1),0,pi*ones(1,NIR),acqParam.FA_true_rads(acqParam.isIR==1), ...
            acqParam.NReadout(acqParam.isIR==1),acqParam.PECentre(acqParam.isIR==1)); %signal of IR components
        
        sigma_signal_noise = SI(2)/acqParam.T1_SNR; % standard deviation of noise
        T1_noise = sigma_signal_noise * randn(size(SI)); % random array of noise, same size as SI
        SI = SI + T1_noise; % add noise to signal
        
        % fit T1 scan intensities to model
        [T1_meas,S0_meas,k_meas,modelFit] = fit_R1(SI,acqParam.isIR,acqParam.isFit,acqParam.TR_s,acqParam.FA_nom_rads,acqParam.TI_s,acqParam.PECentre,acqParam.NReadout,acqParam.NTry);
    case 'Assumed'
        assumed_T1 = varargin{1};
        T1_meas = assumed_T1;
        S0_meas = NaN;
        k_meas = NaN;
        modelFit = NaN;
    otherwise
        error('Error: T1 acquisition method not recognised')        
end