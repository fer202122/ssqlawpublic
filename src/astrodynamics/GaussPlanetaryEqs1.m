%% GaussPlanetaryEqs1: Gauss Planetary Equations with classical orbital element set. 
% Take true anomaly as a fast variable
% IN THE SPACECRAFT-CENTERED HILL coordinate frame
function Bsch = GaussPlanetaryEqs1(GM, OrbElem)
a = OrbElem(1);
e = OrbElem(2);
I = OrbElem(3);
omega = OrbElem(5);
Omega = OrbElem(4);
TA = OrbElem(6); % True anomaly
p = a*(1-e^2);
r = p/(1+e*cos(TA));

h = sqrt((1- e^2) * GM * a); % [km^2/s]

Bsch = h/GM*[2*a/(1-e^2)*e*sin(TA), 2*a/(1-e^2)*(1+e*cos(TA)), 0;... % [s]
			sin(TA), ((e+cos(TA))/(1+e*cos(TA)) + cos(TA)), 0;... % [s/km]
			0, 0, cos(omega+TA)/(1+e*cos(TA));... % [s/km]
			-cos(TA)/e, (2+e*cos(TA))*sin(TA)/(e*(1+e*cos(TA))), - sin(omega+TA)/(tan(I)*(1+e*cos(TA)));... % [s/km]
			0, 0, sin(omega+TA)/(sin(I)*(1+e*cos(TA)));... % [s/km]
			cos(TA)/e, -(2+e*cos(TA))*sin(TA)/(e*(1+e*cos(TA))), 0]; % [s/km] % TA
end