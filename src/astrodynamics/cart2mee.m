function mee = cart2mee(rv, mu)
%CART2MEE Transforms cartesian coordinates to Modified Equinoctial Elements
r = rv(:,1:3);
v = rv(:,4:6);

rdv      = vecDotVec(r,v);
rhat     = makeUnit(r);
rmag     = vecNorm(r);
hvec     = vecCrossVec(r,v);
hmag     = vecNorm(hvec);
hhat     = makeUnit(hvec);
vhat     = (rmag.*v - rdv.*rhat)./hmag;
p        = hmag.*hmag ./ mu;
k        = hhat(:,1)./(1.0 + hhat(:,3));
h        = -hhat(:,2)./(1.0 + hhat(:,3));
kk       = k.*k;
hh       = h.*h;
s2       = 1.0+hh+kk;
tkh      = 2.0.*k.*h;
ecc      = vecCrossVec(v,hvec)./mu - rhat;
fhat(:,1)  = 1.0-kk+hh;
fhat(:,2)  = tkh;
fhat(:,3)  = -2.0.*k;
ghat(:,1)  = tkh;
ghat(:,2)  = 1.0+kk-hh;
ghat(:,3)  = 2.0.*h;
fhat     = fhat./s2;
ghat     = ghat./s2;
f        = vecDotVec(ecc,fhat);
g        = vecDotVec(ecc,ghat);
L        = atan2(rhat(:,2)-vhat(:,1),rhat(:,1)+vhat(:,2));

mee = [p,f,g,h,k,L];

end

