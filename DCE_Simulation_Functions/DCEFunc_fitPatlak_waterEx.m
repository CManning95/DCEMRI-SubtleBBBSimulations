function [PKP, enhModelFit_pct]=DCEFunc_fitPatlak_waterEx(tRes_s,enh_pct,Cp_AIF_mM,Hct,T10_tissue_all_s,T10_blood_s,TR_s,TE_s,FA_all_deg,r1_permMperS,r2s_permMperS,opts)
% Fit concentration curve using model to generate pharmacokinetic
% parameters. At present, only the Patlak model is supported but other
% models may be added on request.
% Output:
% PKP = struct containing pharmacokinetic parameters (vP, PS_perMin, vE); each parameter is a row vector, with each entry
% corresponding to a time series
% enhModelFit_pct = best-fit enhancement curve; each column corresponds to
% a time series.
% Input:
% tRes_s = time resolution in seconds for input data
% enh_pct: array of enhancements; each column corresponds to a time series
% Cp_AIF_mM = column vector giving AIF plasma (not blood) concentration in mM
% Hct = haematocrit
% T10_tissue_all_s: row vector containing pre-contrast tissue T1 values for each time series
% T10_blood_s: blood T1 value
% TR_s, TE_s: acquisition parameters for SPGR/FLASH sequence
% FA_all_deg: row vector of flip angles; each value corresponds to a time series
% r1_permMperS, r2_permMperS: T1 and T2 relaxivity values. Ignores T2'
% effects.
% opts = struct containing options:
%   NIgnore = number of points to exclude in cost function. For Patlak,
%       early data points may be excluded to reduce CBF effects
%   init_vP = initial value for vP
%   init_PS_perMin = initial value for PS

N=size(enh_pct,2); %number of time series
NTime=size(enh_pct,1); %number of time points

enhModelFit_pct=nan(NTime,N); %initialise outputs
PKP.vP=nan(1,N);
PKP.PS_perMin=nan(1,N);

fitParams0=[opts.init_vP opts.init_PS_perMin]; %vector of starting values for fit

T2s0_tissue_s=1; %nominal value (value does not affect enhancement)

for iSeries=1:N %loop through different time series (e.g. voxels or ROIs)
    if sum(isnan(enh_pct(:,iSeries))) > 0; continue; end % skip concentration profiles containing one or more NaNs
    
    %% set T10_tissue and FA, which may be different for each time series
    T10_tissue_s=T10_tissue_all_s(1,iSeries);
    FA_deg=FA_all_deg(1,iSeries);
    
    %% Optimise vP and PS until predicted enhancement matches measured enhancement
    [fitParams,resnorm,residual,exitflag,output]=lsqcurvefit(@fPatlak,fitParams0,[],enh_pct(opts.NIgnore+1:end,iSeries)...
        ,[0 -inf],[1 inf],optimoptions(@lsqcurvefit,'Display','iter','TolFun',1e-7,'TolX',1e-7,'TypicalX',fitParams0));
    
    if exitflag<1
        continue; %return nan if solution not found
    end
    
    enhModelFit_pct(:,iSeries)=[nan(opts.NIgnore,1); fPatlak(fitParams,[])]; %return enhancement curve for best fit parameters
    PKP.vP(1,iSeries)=fitParams(1); %get best fit parameters
    PKP.PS_perMin(1,iSeries)=fitParams(2);
end

%% Function to calculate enhancement for any combination of vP and PS. This is the function that will be fitted to the enhancement curve
    function enhModel_pct=fPatlak(p,x) %function to minimise when using non-linear implementation
        testPKParams.vP=p(1);
        testPKParams.vB=p(1)/(1-Hct);
        testPKParams.PS_perMin=p(2);
        testPKParams.vE=1-testPKParams.vB; %For Patlak, define vE to cover all non-blood space. This means c_e_mM will represent the mean concentration across all extravascular space
        
        %We assume EES and I spaces are in fast water exchange, i.e. acting as a single compartment
        %Practically, we model this by setting [vB vE vI] -> [vB 1-vB 0] (i.e. we use EES to represent the entire extravascular compartment)
        
        %v = volume fractions [vB vE vI]
        v = [testPKParams.vB testPKParams.vE 0 ].';
        p = v; %assume population fractions are same as volume fractions
        
        %Use Patlak model to determine capillary blood and EES concentrations...
        [Ct_mM, IRF, c_cp_mM, c_e_mM] = DCEFunc_PKP2Conc(tRes_s,Cp_AIF_mM,testPKParams,'Patlak');
        %Calculate T10 for EES based on blood and overall tissue values, assuming the FXL...
        T10_ES_s = (1-p(1)) * ((T10_tissue_s*T10_blood_s)/(T10_blood_s-p(1)*T10_tissue_s));
        T10_s = [T10_blood_s T10_ES_s T10_ES_s].'; % vector of compartmental pre-contrast T10
        %Calculate enhancement assuing the SXL...
        enhModel_pct = H2OEx_Conc2Enh_SPGR((1-Hct)*c_cp_mM, c_e_mM, v, p, T10_s, T2s0_tissue_s,TR_s,TE_s,FA_deg,[],[],r1_permMperS,r2s_permMperS,'SXL');
        enhModel_pct = enhModel_pct(opts.NIgnore+1:end);
        PKP.Ct_SXL_mM(:,iSeries) = Ct_mM;
    end

end


