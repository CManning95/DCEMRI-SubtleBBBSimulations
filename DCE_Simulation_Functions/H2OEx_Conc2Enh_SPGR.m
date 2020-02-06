function enhancementPct=H2OEx_Conc2Enh_SPGR...
    (conc_cp_mM,conc_EES_mM,v,p,T10_s,T2s0_s,TR_s,TE_s,FA_deg,kbe_perS,kie_perS,r1_permMperS,r2s_permMperS,exch_model)

% calculate signal enhancement from contrast agent concentration for SPGR
% sequence taking account of water exchange effects
% OUTPUT:
% enhancementPct: array of enhancements; each column = 1 time series
% INPUT:
% conc_cb_mM: array of local capillary concentrations (i.e. local as in mMol of CA per unit volume of capillary blood); each column = 1 time series
% conc_EES_mM: array of local EES concentrations; each column = 1 time series
% v: array [vb (volume fraction of blood) ve (volume fraction of EES)] - used only for calculating total tissue Gd concentration
% p: array of H2O population fractions in the [blood EES intracellular] compartments
% T10_s: array of pre-contrast T1 values in the [blood EES intracellular] compartments
% T2s0_s: pre-contrast T2* of tissue
% TR_s, TE_s: acquisition parameters for SPGR/FLASH sequence
% FA_deg: row vector of flip angles; each value corresponds to a time series
% kbe_perS: blood-EES water exchange rate
% kie_perS: intracellular-EES water exchange rate
% r1_permMperS, r2_permMperS: T1 and T2 relaxivity values. Ignores T2' effects
% exch_model: string specifying water exchange model to use:
%   3S2X_num: numerical calculation of relaxation components for 3S2X model
%   FXL: assume fast exchange, i.e. single component with weighted average R1
%   SXL: assume slow exchange, i.e. 3 components relax independently
%   2S1XA: 2-site analytical model (Dickie et al., NeuroImage 184 (2019):349, based on previous work by Buckley (2008) and others

% calculate intracellular and tissue concentrations
conc_i_mM = zeros(size(conc_cp_mM)); % intracellular local Gd concentration - set to zero array, as Gd contrast agent cannot permeate cells
conc_t_mM = (v(1)*conc_cp_mM) + (v(2)*conc_EES_mM); % overall tissue Gd concentration - set to weighted sum of capillary and EES concentrations

% calculate initial equilibrium magnetisation vector, i.e. [M0(blood) M0(EES) M0(intracellular)]
M0 = p; % set magnetisation vector equal to water population fraction, since any scale factor is irrelevant for calculating enhancement

%seperates T1s into components
T10_b_s = T10_s(1); %blood
T10_EES_s = T10_s(2); %EES
T10_i_s = T10_s(3); % intracellular

%creates arrays
enhancementPct = nan(size(conc_cp_mM));
NTimePoints=size(conc_cp_mM,1);
NSeries=size(conc_cp_mM,2);

for iSeries=1:NSeries %loop through time series (e.g. voxels or ROIs)
    signal_pre = H2OEx_getSPGRSignal(M0,T10_s,T2s0_s,TR_s,TE_s,FA_deg(1,iSeries),kbe_perS,kie_perS,exch_model); %calculate pre-contrast signal
    for iTimePoint=1:NTimePoints %loop through time points within this series
        T2s_s = ( T2s0_s^-1 + r2s_permMperS * conc_t_mM(iTimePoint) ).^-1;... %calculate T2* post-enhancement based on overall tissue Gd concentration
            
        T1_s = [(T10_b_s^-1 + r1_permMperS * conc_cp_mM(iTimePoint) ).^-1; ... %calculate T1 post-enhancement based on local Gd concentration in each compartment
                (T10_EES_s^-1 + r1_permMperS * conc_EES_mM(iTimePoint) ).^-1; ...
                (T10_i_s^-1 + r1_permMperS * conc_i_mM(iTimePoint) ).^-1];
            
        signal_post = H2OEx_getSPGRSignal(M0,T1_s,T2s_s,TR_s,TE_s,FA_deg(1,iSeries),kbe_perS,kie_perS,exch_model); %calculate post-contrast signal based on the modified T1 and T2 values
        enhancementPct(iTimePoint) = 100 * ( (signal_post - signal_pre) / signal_pre ); %calculate signal enhancement
    end
end

        
