function output = coneOpt(psi)
%CONEOPT Computes optimal cone angle that yields maximum ss-acc projected along unit vector with cone angle Psi
output = 1.0/2.0*(psi - asin(sin(psi)/3.0)); % [-pi/2, pi/2]

end