function output = sunSincInc(coe)
% SUNSININC Computes inclination in ECI that yields Sun-Synchronous raan drift
% Constants
const = getConstants();
J2 = const.J2; % [-]
Re = const.Re; % [km]
mu = const.mu; % [kmˆ3/sˆ2]
OmegaDot = deg2rad(360/1*1/const.syear); % RAAN drift to match the apparent motion of the Sun [rad/s]

% SMA, ecc and inclination
sma = coe.sma;
ecc = coe.ecc;

% Orbital angular velocity
n = sqrt(mu/sma^3); % [rad/s]

% Sun-synchronous inclination
A = -3/2*J2*(Re/(sma*(1 - ecc^2)))^2*n; % [rad/s]
output = acos(OmegaDot/A); % [rad]

end