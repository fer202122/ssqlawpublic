function aa = aspectAngle(cart_P, sunHat_P)
%ASPECTANGLE Computes the evolution of the Aspect Angle (AA)
r_P = cart_P(:,1:3);
v_P = cart_P(:,4:6);

% Angular momentum unit vector in ECI frame (i.e., hhat_E)
hhat_P = vecCrossVec(r_P, v_P);
hhat_P = hhat_P./vecNorm(hhat_P);

% Aspect angle at t0 = 0
aa = real(acos(vecDotVec(hhat_P,sunHat_P)));

end