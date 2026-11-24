function output = dJdcone(alpha, phi, params)
% DJDCONE Partial derivative of J = a_s*lambda used to compute coneOpt
C1 = params.C1;
C2 = params.C2;
C3 = params.C3;

% Compute value of partial derivative
k = cos(phi)./sin(phi);
output = k.*sin(alpha).*(3*C1*cos(alpha).^2 + 2*C2*cos(alpha) + C3)...
         - C1*cos(alpha).*(1 - 3*sin(alpha).^2) - C2*cos(2*alpha);

end