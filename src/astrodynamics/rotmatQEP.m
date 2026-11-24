function QEP = rotmatQEP(obl)
arguments
    obl = getConstants().obl;
end

% Rotation around the x-axis
QEP = Rot_body(1, -obl);

end

