function [T1_components_s, w, w_inv, M0prime] =...
    H2OEx_getT1Components_3S2X_num(M0,T1_s,kbe_perS,kie_perS)
% Calculates the T1 components and equilibrium magnetisations assuming a 3-site, 2-exchange model.
% INPUT:
% vectors list components as follows: [ Blood Extracellular Intracellular]
% M0: equilibrium magnetisation for each component (column vector)
% T1_s: T1 for each compartment
% kbe_perS: exchange rate constant from blood to EES
% kie_perS: exchange rate constant from intracellular to EE space
% OUTPUT:
% T1_components_s: T1 for each relaxation component
% w: matrix whose columns are the eigenvectors of A; transforms
% magnetisation from eigenbasis to physical basis: M = w*M'
% w_inv: inverse of w and transforms magnetisation from physical basis to
% relaxation components basis: M'=w_inv*M. i.e. a row of w_inv represents the
% contribution of blood, extraC and intraC magnetisation for a given
% component

%% set equilibrium magnetisations
M0b = M0(1);
M0e = M0(2);
M0i = M0(3);

%% use mass balance to calculate additional rate constants
keb = (kbe_perS*M0b)/M0e;
kei = (kie_perS*M0i)/M0e;

%% set relaxation rates
R1b=1/T1_s(1); R1e=1/T1_s(2); R1i=1/T1_s(3);

%% Set out Bloch equation
% dM/dt = A.M + C
% where M = [ Mb; Me; Mi ];
A = [ -R1b-kbe_perS        keb        0            ;...
         kbe_perS     -R1e-keb-kei   kie_perS      ;...
            0              kei     -R1i-kie_perS   ];

% C = [ M0b*R1b ;...
%       M0e*R1e ;...
%       M0i*R1i ];

%% diagonalise A, yielding the eigenvalues R1_1, R1_2, R1_3
[w_unsorted,A_diag] = eig(A); %get eigenvalues and eigenvectors of A
[T1_components_s,idx] = sort(-diag(1./A_diag)); %sort in order of T1
w = w_unsorted(:,idx);

% rows in w_inv represent the proportion of blood, extra and intra signal for each relaxation component
w_inv = w\eye(3,3);

% calculate M0 in new basis
M0prime = w_inv * M0;

end