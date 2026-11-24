function mee = coe2mee(coe)
a = coe(:,1);
ecc = coe(:,2);
inc = coe(:,3);
Ome = coe(:,4);
ome = coe(:,5);
nu = coe(:,6);

p = a.*(1-ecc.^2);
f = ecc.*cos(ome + Ome);
g = ecc.*sin(ome + Ome);
h = tan(inc./2).*cos(Ome);
k = tan(inc./2).*sin(Ome);
L = Ome + ome + nu;

mee = [p,f,g,h,k,L];

end
