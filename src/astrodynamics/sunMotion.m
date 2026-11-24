function thetas = sunMotion(t, t0jd, const)
arguments
    t
    t0jd
    const = getConstants();
end
% Delta time in Julian Date
dtjd = t/const.sday; % [jd]
% Time in Julian Date
tjd = t0jd + dtjd; % [jd]
% Ephemeris function call
auxdata.const = const;
[~, thetas] = sunMotionEclipticCircular(tjd, auxdata);

end

