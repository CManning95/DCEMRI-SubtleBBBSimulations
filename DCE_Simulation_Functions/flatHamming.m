function w = flatHamming(N,pct);
%Window function that is 1 for the central 100-pct %
%and decays as a Hamming function on either side.
%Centred at N/2 + 1 (assumes even N)

w=nan(1,N);
centre=N/2 + 1;
flatPct=100-pct;
flatPointLeft=round(centre - 0.5*N*(flatPct/100));
if flatPointLeft<1; flatPointLeft=1; end;
flatPointRight=2*centre-flatPointLeft;
if flatPointRight>N; flatPointRight=N; end;

w(flatPointLeft:flatPointRight)=1;

for idx=1:(flatPointLeft-1)
    w(idx) = 0.54 - 0.46 * cos(pi*((idx-1)/(flatPointLeft-1)));
end
w((flatPointRight+1):N)=w((flatPointLeft-1):-1:2);

