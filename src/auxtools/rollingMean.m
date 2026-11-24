function [COEm] = rollingMean(t, COE, mu)
%ROLLINGMEAN Computes rolling mean elements, where the considered window is the current orbital period
%   Detailed explanation goes here
arguments (Input)
    t
    COE
    mu
end

arguments (Output)
    COEm
end

% Store COE history and time
persistent COEs
persistent ts
if t==0
    COEs = [];
    ts = [];
end
COEs = [COEs; COE'];
ts = [ts; t];

% Compute mean elements if there are more than one element in the history
if length(ts) > 1
    % Select time-window within current orbital period
    Tc = 2*pi*sqrt(COE(1)^3/mu);
    idxs = find(ts>t-Tc);
    tw = ts(idxs);
    COEw = COEs(idxs,:);

    % Compute mean SMA and eccentricity
    COEm = zeros(6,1);
    COEm(1:2) = mean(COEw(:,1:2));

    % Compute circular mean for i, Ome, ome [-pi pi]
    COEm(3:5) = atan2(sum(sin(COEw(:,3:5))), sum(cos(COEw(:,3:5))));

    % Obtain true anomaly w.r.t the mean argument of perigee
    u = COE(5) + COE(6); % Osculating argument of latitude ome + true
    COEm(6) = u - COEm(5); % True anomaly w.r.t mean argument of perigee

    % Reduce size of COE history and time to save memory
    COEs = COEw;
    ts = tw;

else
    COEm = COE;

end

end