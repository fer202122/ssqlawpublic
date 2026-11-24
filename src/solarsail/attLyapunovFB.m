%% attLyapunovFB: function description
function [att] = attLyapunovFB(chi, SRP_coefficients)
if numel(chi(:,1))==1
	chi = chi(:);
end
Cx = chi(1,:); Cy = chi(2,:); Cz = chi(3,:);

% k = Cx./sqrt(Cy.^2+Cz.^2);

% alp0 = atan(1./4.*(-3.*k+sqrt(8+9.*k.^2)));
bet = atan2(Cy,Cz);

% FGL: I need to review the signs in these formulas cause I reversed the normal
% -------------------------------------------------------------------------
%{
C1 = SRP_coefficients(1);
C2 = SRP_coefficients(2);
C3 = SRP_coefficients(3);

ca0 = cos(alp0);
sa0 = sin(alp0);

F 		= k.*(3.*C1.*ca0.^2 + 2.*C2.*ca0 + C3).*sa0 - C1.*ca0.*(1-3.*sa0.^2) - C2.*(cos(2.*alp0));
dFda 	= k.*(3.*C1.*(ca0.^3 - 2.*ca0.*sa0.^2) + 2.*C2.*cos(2.*alp0) + C3.*ca0) - C1.*sa0.*(2 - 9.*ca0.^2) + 2.*C2.*sin(2.*alp0);
% F 		= -k.*(3.*C1.*ca0.^2 + 2.*C2.*ca0 + C3).*sa0 + C1.*ca0.*(1-3.*sa0.^2) + C2.*(cos(2.*alp0));
% dFda 	= -k.*(3.*C1.*(ca0.^3 - 2.*ca0.*sa0.^2) + 2.*C2.*cos(2.*alp0) + C3.*ca0) + C1.*sa0.*(2 - 9.*ca0.^2) -2.*C2.*sin(2.*alp0);

alp = alp0 - F./dFda;
%}
% -------------------------------------------------------------------------
% Primer cone (phi) and clock (delta) angles
% Cone. We use the more numerically stable version defined by Oskar
primer_S = chi;
primer_S = makeUnit(primer_S);
normdiff = vecNorm(primer_S - [1; 0; 0]);
normsum = vecNorm(primer_S + [1; 0; 0]);
phi = 2*atan2(normdiff, normsum); % [0 pi]

% Critical primer cone angle (phiCrit) for current sail parameters (from gpops2cases/src)
params.C1  = SRP_coefficients(1);
params.C2  = SRP_coefficients(2);
params.C3  = SRP_coefficients(3);
[~, phiCrit] = coneCrit(params);

% Decide regime based on critical primer cone angle
if phi==0                           % Trivial solution
    alpha = 0.0;
elseif (0 < phi) && (phi < phiCrit) % Oguri's approximation regime
    alpha = coneOptOguri(phi, params);
elseif phiCrit < phi                % Deactivate sail, it is not possible to generate positive thrust along primer
    alpha = pi/2;
end

% Output
att = [alpha; bet];

end
