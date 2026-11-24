function [alphaCrit, phiCrit] = coneCrit(params)
%CONECRIT Computes the critical cone angle and the critical cone angle
% of the primer vector for which the sail has to be deactivated
if ~isstruct(params)
    par.C1  = params(1);
    par.C2  = params(2);
    par.C3  = params(3);
else
    par = params;
end

% Optical parameters
C1 = par.C1;
C2 = par.C2;
C3 = par.C3;

% Rename to match Mengali and Quarta parameters
b1 = C3;
b2 = C1;
b3 = C2;

% Compute critical angles
alphaCrit = (-b1*b3 - 2*b2*b3 + sqrt(b1^2*b3^2 - 4*b1*b3^2*b2 + 8*b1^2*b2^2 + 4*b2^3*b1))/(4*b1*b2 + 2*b2^2);
alphaCrit = acos(alphaCrit);

phiCrit = primerConeOpt(alphaCrit, par);

end