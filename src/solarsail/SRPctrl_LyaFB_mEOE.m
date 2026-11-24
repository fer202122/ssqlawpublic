function [att, auxout] = SRPctrl_LyaFB_mEOE(t, mEOE, mEOEtarg, f_maxAlongPrimer, prms, rkFlag)
% Computes the sail cone and clock angles (as defined in Oguri's paper) to maximize the solar-sail
% acceleration along the primer vector derived from the gradient of the Q-law

% Transform MEE to COE
COE = mee2coe(mEOE')';
COEtarg = mee2coe(mEOEtarg')';

% Compute mean orbital elements
if prms.meanf
    if rkFlag % If we're inside the RK integrator, the mean filter should not include the intermediate steps
        COE = rollingMeanRk(t, COE, prms.mu);
    else
        COE = rollingMean(t, COE, prms.mu);
    end
end

% Compute the gradient of the Q function via finite differences
[Q, dQdx, tau, dtaudx, K, OrbDist, oenndot, Lopt, oenndotAn] = QgainConstGradDeep(t, COE, COEtarg, prms);

% Compute primer in Local frame (L)
Bmat_L = GaussEqClippedCoe(prms.mu, COE); % [s] for SMA and [s/km] for the rest
primer_L = - Bmat_L'*dQdx'; % [s^3/km]

% Compute the sail's attitude that maximizes the controlled acceleration along the primer
att = f_maxAlongPrimer(t, mEOE, primer_L);

% Auxiliary output
auxout.primer_L = primer_L;
auxout.Q = Q;
auxout.dQdx = dQdx;
auxout.Bmat_L = Bmat_L;
auxout.K = K;
auxout.tau = tau;
auxout.OrbDist = OrbDist;
auxout.oenndot = oenndot;
auxout.Lopt = Lopt;
auxout.dtaudx = dtaudx;
auxout.oenndotAn = oenndotAn;

end
