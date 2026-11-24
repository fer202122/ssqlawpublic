function [sunHat_P, thetas] = sunMotionEclipticCircular(tjd, auxdata)
%SUNMOTIONECLIPTICIRCULAR Computes the Sun's unit vector in the Ecliptic Frame
% assuming circular motion for the Sun
const = auxdata.const; % AdiGator does not like getConstants()

N = length(tjd);

t0jd = const.autEquiJd;
omegas = 2*pi/const.syear*const.sday; % [rad/day]
thetas = omegas*(tjd - t0jd);

sunHat_P = [cos(thetas), sin(thetas), zeros(N,1)];

end