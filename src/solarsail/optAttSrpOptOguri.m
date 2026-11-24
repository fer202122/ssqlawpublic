function att = optAttSrpOptOguri(t, mEOE, primer_L, Cprms, f_QLS)
%OPTATTSRPOPTOGURI Optimizes the attitude of the sail that maximizes the projection of the SRP acceleration along the primer vector using the optical model and Oguri's approach

% Rotate primer from Local frame (L) to the Sun frame (S)
coe = mee2coe(mEOE')';
QSL = f_QLS(t, coe)';
primer_S = QSL*primer_L;

% Primer cone (phi) and clock (delta) angles
% Cone. We use the more numerically stable version defined by Oskar
primer_S = makeUnit(primer_S);
normdiff = vecNorm(primer_S - [1; 0; 0]);
normsum = vecNorm(primer_S + [1; 0; 0]);
phi = 2*atan2(normdiff, normsum); % [0 pi]
% Clock
delta = atan2(primer_S(2), primer_S(3)); % [-pi, pi]

% Critical primer cone angle (phiCrit) for current sail parameters (from gpops2cases/src)
params.C1  = Cprms(1);
params.C2  = Cprms(2);
params.C3  = Cprms(3);

% Oguri C1
alpha = coneOptOguri(phi, params);
alpha = median([alpha; [0; pi/2]]);

% Output
att = [alpha; delta];

end
