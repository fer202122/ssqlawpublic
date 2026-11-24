function r = mee2r(mee)
%mee2r Transforms MEE state to position vector magnitude

% MEE
p = mee(:,1);
f = mee(:,2);
g = mee(:,3);
L = mee(:,6);

% Radius vector magnitude
q = 1+f.*cos(L)+g.*sin(L);
r = p./q;

end

