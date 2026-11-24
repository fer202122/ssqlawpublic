function QPS = rotmatQPS(t, t0jd)
% Sun's current position in the Ecliptic frame:
thetas = sunMotion(t, t0jd);

% Minus rotation around z-axis
QPS = Rot_body(3,-1*thetas);

end

