function mEOE_E = mEOEsP2E(mEOE_P, mu)
%MEOEE2P Transforms Modified Equinoctial Elements from Ecliptic (P) to Earth (E) Frame

% Transform mEOE to cartesian
auxdata.mu = mu;
cart_P = mee2cart(mEOE_P, auxdata);


% Rotate cartesian state from Ecliptic (P) to Earth (E) frame
QEP = rotmatQEP;
r_P = cart_P(:,1:3);
v_P = cart_P(:,4:6);
r_E = matDotVec(QEP, r_P);
v_E = matDotVec(QEP, v_P);
cart_E = [r_E,v_E];

% Transform to mEOE
mEOE_E = cart2mee(cart_E, auxdata.mu);

end