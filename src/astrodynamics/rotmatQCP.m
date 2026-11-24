function QCP = rotmatQCP(coe_P)
%ROTMATQCP Rotation matrix from Ecliptic to Cranked
raan = coe_P.raan;
inc = coe_P.inc;
N = length(raan);

QCP = cat(3,[           cos(raan),            sin(raan), zeros(N,1)],...
            [-cos(inc).*sin(raan),  cos(inc).*cos(raan),   sin(inc)],...
            [ sin(inc).*sin(raan), -sin(inc).*cos(raan),   cos(inc)]);
QCP = permute(QCP, [3,2,1]);

end