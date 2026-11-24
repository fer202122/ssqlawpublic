function output = dJdcone2(alpha, phi, params)
% DJDCONE Second partial derivative of J = a_s*lambda used to compute coneOpt
C1 = params.C1;
C2 = params.C2;
C3 = params.C3;

% Compute value of partial derivative
k = cos(phi)./sin(phi);
output = k.*(3*C1*(cos(alpha).^3 - 2*cos(alpha).*sin(alpha).^2) + 2*C2*cos(2*alpha) + C3*cos(alpha))...
         - C1*sin(alpha).*(2 - 9*cos(alpha).^2) + 2*C2*sin(2*alpha);

end