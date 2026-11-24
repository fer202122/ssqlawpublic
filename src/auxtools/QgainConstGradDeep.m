function [Q, dQdx, tau, dtaudx, K, OrbDist, oenndot, Lopt, oenndotAn] = QgainConstGradDeep(t, COE, COEtarg, prm)
%QgainNumGradDeep Compute dQdx by finite differentiation but at a deeper level (no Q-function required)
% Check for negative eccentricity
if COE(2) < 0
    fprintf('Negative e at time %1.3f [days] \n', t/86400);
    COE(2) = 0;
end

% oenndot values (constant)
output = COEnnConst(t, COE, COEtarg, prm);
oenndot = output.oenndot;
oenndotAn = output.oenndotAn;
Lopt = output.Lopt;

% oenndot values for comparison:
% % oenndot numerical values for comparison
% output = COEnnNum(t, COE, COEtarg, prm);
% oenndotAn = output.oenndot;
% oenndot analytical values for comparison
output = COEnnAna(t, COE, COEtarg, prm);
oenndotAn = output.oenndot;

% Safeguard for oenndot values close to zero or positive
oenndot = safeguard(oenndot, COEtarg, prm);

% Evaluate distance function
OrbDist = coeDistance(COE, COEtarg);

% tau vector in sma and eccentricity
tau = -1*abs(OrbDist)./oenndot;

% tau partials matrix
sig = cmptSign(coeDistance(COE, COEtarg));
dtaudx = zeros(2,2);
dtaudx(:,1) = (-1)/oenndot(1)*[ sig(1) ;...
                                    0  ];
dtaudx(:,2) = (-1)/oenndot(2)*[     0  ;...
                                sig(2) ];

% Interpolate weights if time-dependet
if isstruct(prm.W)
    auxdata.control = prm.W;
    weights = uInterp(t, [], auxdata);
else
    weights = prm.W;
end

% Assemble Q
K = diag(weights);
Q = tau'*K*tau;

% Assemble Q-gradient
dQdx = 2*dtaudx*K*tau;
dQdx = [dQdx', zeros(1,4)]; % output (we need a row vector)

end