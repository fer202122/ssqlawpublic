function phiOpt = primerConeOpt(alpha, params)
%PRIMERCONEOPT Computes the optimal primer cone angle for a given normal 
% cone angle (optical model)

% Optical parameters
C1 = params.C1;
C2 = params.C2;
C3 = params.C3;

% Compute optimal primer cone angle
phiOpt = sin(alpha).*(3*C1*cos(alpha).^2 + 2*C2*cos(alpha) + C3)./...
                    (cos(alpha).^2.*(C1*cos(alpha) + C2) - sin(alpha).^2.*(2*C1*cos(alpha) + C2));
phiOpt = atan(phiOpt);

% Wrap to [0 pi]
phiOpt(phiOpt<0) = phiOpt(phiOpt<0) + pi;

end
