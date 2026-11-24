function T = orbPeriod(sma, mu)
%ORBPERIOD Computes the orbital period
T = 2*pi*sqrt(sma.^3/mu);

end