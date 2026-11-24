function signOrbDist = cmptSign(orbDist)
%CMPTSIGN Computes the sign of the orbital distance
signOrbDist = sign(orbDist);

% If the orbital distance is zero, we assign a positive sign
signOrbDist(signOrbDist==0) = 1;

end