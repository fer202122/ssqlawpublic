function mEOE_P = mEOEsE2P(mEOE_E, mu)
%MEOEE2P Transforms Modified Equinoctial Elements from Earth (E) to Ecliptic (P) Frame

% Transform mEOE to cartesian
auxdata.mu = mu;
cart_E = mee2cart(mEOE_E, auxdata);


% Rotate cartesian state from Earth (E) to Ecliptic (P) frame
QPE = rotmatQEP';
r_E = cart_E(:,1:3);
v_E = cart_E(:,4:6);
r_P = matDotVec(QPE, r_E);
v_P = matDotVec(QPE, v_E);
cart_P = [r_P,v_P];

% Transform to mEOE
mEOE_P = cart2mee(cart_P, auxdata.mu);

end
