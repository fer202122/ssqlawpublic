function aa = wrapAa(aa)
%WRAPAA Wrap the aa to the [0,pi/2] interval
aa(aa>pi/2) = pi - aa(aa>pi/2);

end

