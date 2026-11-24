function u = uInterp(t, x, auxdata)
% Interpolate control
    if auxdata.control.dim == 1
        u = ppval(auxdata.control, t);
    else
        u = ppval(auxdata.control, t)';
    end

end
