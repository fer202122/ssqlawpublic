function alphaOpt = coneOptOp(phi, params)
%CONEOPTOP Optimal cone angle (optical) for a given primer cone angle
% Critical primer cone angle
[~, phiCrit] = coneCrit(params);

% Initialize as pi/2
alphaOpt = pi/2*ones(numel(phi),1);

% For those values smaller than the critical, compute optimum cone angle
alphaOpt(phi < phiCrit) = coneOptOguri(phi(phi < phiCrit), params);
% alphaOpt(phi < phiCrit) = coneOptNum(phi(phi < phiCrit), params);

end