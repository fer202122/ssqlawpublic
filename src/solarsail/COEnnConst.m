%% QlawOEnnAna: function description
function output = COEnnConst(t, COE, COEtarg, prm)
% COEnnConst Q-law orbital element double minimization (i.e., nn) constant values
% Initialize output arrays
oenndot = zeros(2,1);
oenndotAn = zeros(2,1);
TA_opt  = zeros(2,1);

% Constants
a0   = 2*prm.g0;        % Characteristic acceleration
sma0 = prm.COE0(1);
% inc0 = prm.COE0(3);
h0 = sqrt(sma0*prm.mu); % Initial angular momentum magnitude

%% Semi-major axis oenn
% % Primer vector has a cone angle phi1 wrt the Sun direction
% % The optimal normal orientation that maximizes the projection of the 
% % solar-sail acceleration along the primer vector is computed
% phi1 = inc0;
% % alphaOpt = coneOpt(phi1);
% % oenndot(1) = (-1/h0*2*sma0^2)*(a0*cos(alphaOpt)^2)*cos(phi1 - alphaOpt);

if all(prm.Cprms==[2 0 0])
    % Normal vector and primer vector are parallel (ideal model)
    oenndot(1) = (-1/h0*2*sma0^2)*a0;

else
    % Normal vector and primer vector are parallel (optical model)
    acc_opt 	= a0/2 * work_aSRP_hat([0; 0], prm.Cprms);
    oenndot(1) = (-1/h0*2*sma0^2)*acc_opt(1);

end

% Analytical output is the same as the constant output
oenndotAn(1) = oenndot(1);


%% Eccentricity oenn
% % Primer vector opposed to the primer vector with cone angle phi1
% phi2 = pi - inc0;
% alphaOpt = coneOpt(phi2);
% oenndot(2) = (-1/h0*2*sma0)*(a0*cos(alphaOpt)^2)*cos(phi2 - alphaOpt);

if all(prm.Cprms==[2 0 0])
    % Normal vector and primer vector are perpendicular (ideal model)
    oenndot(2) = (-1/h0*2*sma0)*(a0*cos(prm.const.coneOpt)^2)*sin(prm.const.coneOpt);

else
    % Normal vector and primer vector are perpendicular (optical model)
    p_oe = (1/h0*2*sma0)*[0 0 1];
    att_opt 	= attLyapunovFB(p_oe, prm.Cprms);
    acc_opt 	= a0/2 * work_aSRP_hat(att_opt, prm.Cprms);
    oenndot(2) = -p_oe*acc_opt;

end

%% Analytical output is the same as the constant output
oenndotAn(2) = oenndot(2);

%% output
output.oenndot = oenndot;
output.oenndotAn = oenndotAn; 
output.Lopt = TA_opt;

% %% Test
% annn = -a0/2 * (1/h0*2*sma0^2) * sum(prm.Cprms);
% cone = att_opt(1); C1 = prm.Cprms(1); C2 = prm.Cprms(2);
% ennn = -a0/2 * (1/h0*2*sma0  )*cos(cone)*sin(cone)*(C1*cos(cone) + C2);

end