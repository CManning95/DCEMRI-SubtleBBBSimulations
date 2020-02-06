function Stotal=H2OEx_getSPGRSignal(M0,T1_s,T2s_s,TR_s,TE_s,FA_deg,kbe_perS,kie_perS,exch_model)
%calculate signal for SPGR sequence for a 3-compartment tissue
%OUTPUT:
%Stotal: total signal
%INPUT:
%M0: column vector of equilibrium magnetisations (blood, extra, intra)
%T1_s: T1 for each compartment (blood, extra, intra)
%T2s_s: T2* relaxation time (assumes a single value for whole tissue)
%TR_s, TE_s, FA_deg: acquisition parameters
%kbe_perS: exchange rate constant for transfer of water from blood to EES
%kie_perS: exchange rate constant for transfer of water from intracellular to EE space
%exch_model: string specifying exchange model. Permitted values:
%   3S2X_num: numerical calculation of relaxation components for 3S2X model
%   FXL: assume fast exchange, i.e. tissue relaxes monoexponentially with weighted average R1
%   SXL: assume slow exchange, i.e. the 3 tissue components relax independently
%   2S1XA: 2-site analytical model (Dickie et al., NeuroImage 184 (2019):349, based on previous work by Buckley (2008) and others

%% calculate signal depending on exchange model selected
switch exch_model
    case '3S2X_num' %use numerical full 3S2X model to calculate T1 components, which are the eigenvalues of relaxation/exchange matrix
        [T1_components_s, w, w_inv, M0prime] = H2OEx_getT1Components_3S2X_num(M0,T1_s,kbe_perS,kie_perS);
        
        Sprime = M0prime .* ... %calculate the 3 signal components using the standard SPGR signal equation
            (1-exp(-TR_s./T1_components_s)) ./ (1-cos(2*pi*(FA_deg/360))*exp(-TR_s./T1_components_s)) * ...
            sin(2*pi*(FA_deg/360)) * ...
            exp(-TE_s/T2s_s);
        
        S = w * Sprime; %transform signal back into physical (blood, extra, intra) basis
        Stotal = sum(S,1); %calculate total signal
        
    case 'FXL' %signal behaves as if in a single compartment
        M0total = sum(M0,1); %calculate total magnetisation
        R1eff_perS = sum( (M0 .* (1./T1_s)) , 1)/M0total ; %calculate weighted average R1
        T1eff_s = 1/R1eff_perS; % tissue T1
       
        Stotal =  M0total * ... %calculate single signal component
            (((1-exp(-TR_s./T1eff_s))*sin(2*pi*(FA_deg/360))) / (1-exp(-TR_s/T1eff_s)*cos(2*pi*(FA_deg/360))) ) ...
            * exp(-TE_s./T2s_s) ;
        
    case 'SXL' %calculate separate signal for each compartment, then sum over compartments
        S = M0 .* ... 
            (((1-exp(-TR_s./T1_s))*sin(2*pi*(FA_deg/360))) ./ (1-exp(-TR_s./T1_s)*cos(2*pi*(FA_deg/360))) ) ...
            .* exp(-TE_s./T2s_s);
        Stotal = sum(S,1); %calculate total signal
        
    case '2S1XA'
        Me = M0(2)+M0(3); %equilibrium magnetisation for a single extravascular compartment
        Mi = M0(1); %equilibrium magnetisation for intravascular component
        R1e_perS = ( M0(2)*(1/T1_s(2)) + M0(3)*(1/T1_s(3)) ) / Me; %extravascular R1 is weighted average of extra- and intra-cellular values
        T1e_s = 1/R1e_perS;
        T1b_s = T1_s(1);
        
        pb = Mi/(Mi+Me); %blood population fraction
        tb_s = 1 / kbe_perS; %blood residence time
        
        aS = 0.5 - 0.5 * (... %weighting of short T1 component
            ( (1/T1e_s - 1/T1b_s)*(2*pb-1) + pb/((1-pb)*tb_s) + (1/tb_s) ) / ...
            ( (1/T1e_s - 1/T1b_s + pb/((1-pb)*tb_s) - 1/tb_s )^2 + (4*pb)/((1-pb)*tb_s^2) )^0.5 ...
            );
        
        T1S_s = 2 / ( ... %short T1 component
            (1/T1e_s + 1/T1b_s + pb/((1-pb)*tb_s) + 1/tb_s ) + ...
            ( (1/T1e_s - 1/T1b_s + pb/((1-pb)*tb_s) - 1/tb_s)^2 + (4*pb)/((1-pb)*tb_s^2) )^0.5 ...
            );
        
        T1L_s = 2 / ( ... %long T1 component
            (1/T1e_s + 1/T1b_s + pb/((1-pb)*tb_s) + 1/tb_s ) - ...
            ( (1/T1e_s - 1/T1b_s + pb/((1-pb)*tb_s) - 1/tb_s)^2 + (4*pb)/((1-pb)*tb_s^2) )^0.5 ...
            );
        
        Stotal = (Mi+Me) * (... %total signal
            aS * ( (sin(2*pi*(FA_deg/360))*(1 - exp(-TR_s/T1S_s))) / (1 - cos(2*pi*(FA_deg/360))*exp(-TR_s/T1S_s)) ) ...
            + (1-aS) * ( (sin(2*pi*(FA_deg/360))*(1 - exp(-TR_s/T1L_s))) / (1 - cos(2*pi*(FA_deg/360))*exp(-TR_s/T1L_s)) ) ...
            )...
            *exp(-TE_s./T2s_s);
        
    otherwise
        error('Exchange model not recognised.');
        
end

end
