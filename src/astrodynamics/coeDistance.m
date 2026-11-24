function dist = coeDistance(COE, COEtarg)
%COEDISTANCE Distance function to judge the distance between the current Orbital Elements and the target Orbital Elements
dist = COE(1:2) - COEtarg(1:2);

end