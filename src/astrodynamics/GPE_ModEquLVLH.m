function Bsch = GPE_ModEquLVLH(GM, OrbElem)
mu = GM;
p = OrbElem(1);
f = OrbElem(2);
g = OrbElem(3);
h = OrbElem(4);
k = OrbElem(5);
L = OrbElem(6); % True longitude

% Auxiliary variables
q = 1+f*cos(L)+g*sin(L);
chi = sqrt(h*h+k*k);
s2 = 1+chi*chi;

% Matrix in R,T,N frame
Bsch = [                 0,                 (2*p/q)*sqrt(p/mu),                                     0;...
         sqrt(p/mu)*sin(L), +sqrt(p/mu)*(1/q)*((q+1)*cos(L)+f), -sqrt(p/mu)*(g/q)*(h*sin(L)-k*cos(L));...
        -sqrt(p/mu)*cos(L), +sqrt(p/mu)*(1/q)*((q+1)*sin(L)+g), +sqrt(p/mu)*(f/q)*(h*sin(L)-k*cos(L));...
                         0,                                  0,  sqrt(p/mu)*(s2*cos(L)/(2*q))        ;...
                         0,                                  0,  sqrt(p/mu)*(s2*sin(L)/(2*q))        ;...
                         0,                                  0,  sqrt(p/mu)*(1/q)*(h*sin(L)-k*cos(L))];

end