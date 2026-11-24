function cart = mee2cart(mee, auxdata)
%mee2cart Transforms MEE state to Cartesian state
mu = auxdata.mu;

% MEE
p = mee(:,1);
f = mee(:,2);
g = mee(:,3);
h = mee(:,4);
k = mee(:,5);
L = mee(:,6);

% Auxiliary variables
q = 1+f.*cos(L)+g.*sin(L);
r = p./q;
alpha2 = h.*h-k.*k;
chi = sqrt(h.*h+k.*k);
s2 = 1+chi.*chi;

% Position and velocity vectors
rX = (r./s2).*(cos(L)+alpha2.*cos(L)+2*h.*k.*sin(L));
rY = (r./s2).*(sin(L)-alpha2.*sin(L)+2*h.*k.*cos(L));
rZ = (2*r./s2).*(h.*sin(L)-k.*cos(L));

vX = -(1./s2).*sqrt(mu./p).*(sin(L)+alpha2.*sin(L)-2*h.*k.*cos(L)+g-2*f.*h.*k+alpha2.*g);
vY = -(1./s2).*sqrt(mu./p).*(-cos(L)+alpha2.*cos(L)+2*h.*k.*sin(L)-f+2*g.*h.*k+alpha2.*f);
vZ = (2./s2).*sqrt(mu./p).*(h.*cos(L)+k.*sin(L)+f.*h+g.*k);

% Output
cart = [rX rY rZ vX vY vZ];

end

