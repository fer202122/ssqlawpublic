function coe = mee2coe(mee)
p = mee(:,1);
f = mee(:,2);
g = mee(:,3);
h = mee(:,4);
k = mee(:,5);
L = mee(:,6);

a = p./(1-f.*f-g.*g);
e = sqrt(f.*f+g.*g);
i = atan2(2.*sqrt(h.*h+k.*k),1-h.*h-k.*k);
ome = atan2(g.*h-f.*k,f.*h+g.*k);
Ome = atan2(k,h);
nu = L-Ome-ome;

coe = [a,e,i,Ome,ome,nu];

end
