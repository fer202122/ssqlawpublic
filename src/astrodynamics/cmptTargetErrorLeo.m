function targetError = cmptTargetErrorLeo(L, mEOE, COE0, COEf, mu, relconv)
%CMPTTARGETERRORLEO Error relative to mean J2 oscillation for LEO satellites
% Transform MEE to COE
t = mEOE(7);
COE = mee2coe(mEOE')';

% Compute mean orbital elements
COE = rollingMean(t, COE, mu);

% SMA error
a = COE(1);
a0 = COE0(1);
af = COEf(1);
da = af - a0;
if abs(af)>0 % Relative error w.r.t delta SMA
    aError = abs(a - af)/abs(da);
else % Absolute error
    aError = abs(a - af);
end

% Eccentricity error
e = COE(2);
e0 = COE0(2);
ef = COEf(2);
de = ef - e0;
if ef>0
    eError = abs(e - ef)/abs(de);
else
    eError = abs(e - ef);
end

% Compute maximum relative error and apply the scaling of the relative convergence
targetError =  max([aError; eError].*relconv');
% fprintf('aError, eError, targetError: %1.2e, %1.2e, %1.2e \n', aError, eError, targetError)

end