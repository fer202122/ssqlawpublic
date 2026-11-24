function [ic, betac] = incCrit(params)
%INCCRIT Computes critical inclination below which ecc cannot be controlled no matter the aspect angle
% Prepare input parameters
if ~isstruct(params)
    par.C1  = params(1);
    par.C2  = params(2);
    par.C3  = params(3);
else
    par = params;
end

[~, phiCrit] = coneCrit(par);
cosic = -3*sqrt(2)/4*cos(phiCrit);
ic = acos(cosic);
betac = pi/2 - ic;

end

