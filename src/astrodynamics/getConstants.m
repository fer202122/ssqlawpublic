function output = getConstants()
%GETCONSTANTS Define constant values
% Trial planet:
% % Mars
% planet = 'Mars';
% Re = 3389.5; % Mean radius [km]
% mu = 0.042828e6; % Gravitational parameter [kmˆ3/sˆ2]
% syear = 5.936101426190121e+07; % Sidereal year [s]
% % Jupiter
% planet = 'Jupiter';
% Re = 69911.0; % Mean radius [km]
% mu = 126.687e6; % Gravitational parameter [kmˆ3/sˆ2]
% syear = 3.746235920959806e+08; % Sidereal year [s]
% Earth
planet = 'Earth';
Re = 6378; % Mean radius [km]
mu = 3.986004418e5; % Gravitational parameter [kmˆ3/sˆ2]
syear = 365.256*24*3600; % Earth's sidereal year [s]

sf1AU = 1367.0;        % Solar flux at one Astronomical Unit (AU) [W/m^2]
c     = 299792458.0;   % Speed of light in vacuum [m/s]
srp1AU = sf1AU/c*1e3;  % Solar Radiation Pressure at 1 AU [kg/(km*s^2)]
obl   = deg2rad(23.44);% Obliquity of the Ecliptic [rad]
J2    = 1.08262668e-3; % [-]
rGEO = ((23.9345*3600/(2*pi))^2*mu)^(1/3); % Radius of geostationary orbit [km]
rMEO = 15000 + Re; % User defined radius of MEO orbit [km]
rLEO = 800 + Re; % User defined radius of LEO orbit [km]
autEquiJd = juliandate(datetime(2022,9,23,0,0,0)); % Autum equinox [datetime]
sday = 86400.0; % Seconds in a day [s]

% Planets of the solar system
planets = getPlanets();

% Solar sail
coneOpt = atan(1/sqrt(2)); % Optimal cone angle to maximize acceleration perpendicular to Sun-sail line [rad]
critcAa = pi/2 - 2*atan(1/3*(sqrt(10) - 1.0)); % First critical inclination: yaw discontinuities at hinge  [rad]
critcInc2 = deg2rad(60.0);  % Second critical inclination for which the optimization changes [rad]
% critcInner = deg2rad(65.8605); % Critical inner angle: inner angle of the hinge for the contour at the critical inclination
critcInner = atan(sin(coneOpt)/(cos(coneOpt)*sin(critcAa))); % Critical inner angle: inner angle of the hinge for the contour at the critical inclination

% Geodetic coordinates of W12 (from goolge maps)
w12coord.latitude = 51.990832; % [deg]
w12coord.longitude = 4.371615; % [deg]

% Epochs associated with minimum and maximum solar activity
epochMin = datetime(2019,12,1,0,0,0); % Minimum (December 2019)
epochMax = datetime(2001,12,1,0,0,0); % Maximum (December 2001)

% eta-Psi fitted parameters (i.e., eta = n*Psi^m) for AA = [0, 40, 90]:
% m parameter
m0 = 1.00074654e+00; % [-]
m40 = 1.00062885e+00; % [-]
m90 = 1.00007955e+00; % [-]
% n parameter
n0 = 4.87404000e+00; % [-]
n40 = 2.19795377e+00; % [-]
n90 = 2.39135874e-01; % [-]

% Mapping function (chi) 8th degree polynomial
chiPoly = getChiPoly();

% Maximum Psi
maxPsi = 1e-2; % [-]

% dnu-dgam fitted parameter
% s = 2*const.m0 - 1; % 1.001493080000000: Considering variation of eta_ref with Psi
% s = 2*const.m0 - 5/2; % -0.498506920000000: Considering also the variation of T0 with r
% s = 0.233224635434725; % I don't remember
% s = 0.235929314360624; % From 'multrev_mixedgeo' folder. Fitted to raan=90 & inc==90 curve (best)
% s = 0.474697455806216; % From 'multrev_20241104-175415' folder. Fitted to raan=90 & inc==90 curve (best)
% s = 4.960575917391934e-01; % From 'test/test6' folder. Fitted to raan=90 & inc==90 curve (best)
s = 4.951700036111985e-01; % From 'test/test7' folder. Fitted to raan=90 & inc==90 curve (best)

% Flown sail parameters
sails = getSailsParams(srp1AU);

% drho(dnu) (eta) model polynomial coefficients
drhoPoly = getDrhoPoly();

% output
output.planet = planet;
output.Re = Re;
output.mu = mu;
output.planets = planets;
output.sf1AU = sf1AU;
output.c = c;
output.srp1AU = srp1AU;
output.obl = obl;
output.J2 = J2;
output.rGEO = rGEO;
output.rMEO = rMEO;
output.rLEO = rLEO;
output.syear = syear;
output.autEquiJd = autEquiJd;
output.sday = sday;
output.coneOpt = coneOpt;
output.critcAa = critcAa;
output.critcInc2 = critcInc2;
output.critcInner = critcInner;
output.w12coord  = w12coord;
output.epochMin = epochMin;
output.epochMax = epochMax;
output.m0 = m0;
output.m40 = m40;
output.m90 = m90;
output.n0 = n0;
output.n40 = n40;
output.n90 = n90;
output.maxPsi = maxPsi;
output.s = s;
output.chiPoly = chiPoly;
output.sails = sails;
output.drhoPoly = drhoPoly;

end

function planets = getPlanets()
%GETPLANETS Returns physical constants of various planets of the solar system
% All data was taken from NASA's fact sheets (e.g., https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html (Accesed: 8th of February 2024))
% The Sun
muSun = 132712e6; % Gravitational parameter of the Sun [kmˆ3/sˆ2]
% Mercury
me.R = 2439.7; % Mean radius [km]
me.mu = 0.022032e6; % Gravitational parameter [kmˆ3/sˆ2]
me.sma = 57.909e6; % Semi-major axis [km]
me.soi = soi(me.sma, me.mu, muSun); % Sphere of influence considering the Sun as primary [km]
me.P = orbPeriod(me.sma, muSun); % Orbital period [s]
% Venus
v.R = 6051.8; % Mean radius [km]
v.mu = 0.32486e6; % Gravitational parameter [kmˆ3/sˆ2]
v.sma = 108.210e6; % Semi-major axis [km]
v.soi = soi(v.sma, v.mu, muSun); % Sphere of influence considering the Sun as primary [km]
v.P = orbPeriod(v.sma, muSun); % Orbital period [s]
% Earth
e.R = 6378; % Mean radius [km]
e.mu = 3.986004418e5; % Gravitational parameter [kmˆ3/sˆ2]
e.sma = 149.598e6; % Semi-major axis [km]
e.soi = soi(e.sma, e.mu, muSun); % Sphere of influence considering the Sun as primary [km]
e.P = 365.256*24*3600; % Earth's sidereal year (to be consistent with const.syear) [s]
% e.P = orbPeriod(e.sma, muSun); % Orbital period [s]
% Mars
ma.R = 3389.5; % Mean radius [km]
ma.mu = 0.042828e6; % Gravitational parameter [kmˆ3/sˆ2]
ma.sma = 227.956e6; % Semi-major axis [km]
ma.soi = soi(ma.sma, ma.mu, muSun); % % Sphere of influence considering the Sun as primary [km]
ma.P = orbPeriod(ma.sma, muSun); % Orbital period [s]
% Jupiter
j.R = 69911.0; % Mean radius [km]
j.mu = 126.687e6; % Gravitational parameter [kmˆ3/sˆ2]
j.sma = 778.479e6; % Semi-major axis [km]
j.soi = soi(j.sma, j.mu, muSun); % Sphere of influence considering the Sun as primary [km]
j.P = orbPeriod(j.sma, muSun); % Orbital period [s]
% Saturn
s.R = 58232.0; % Mean radius [km]
s.mu = 37.931e6; % Gravitational parameter [kmˆ3/sˆ2]
s.sma = 1432.041e6; % Semi-major axis [km]
s.soi = soi(s.sma, s.mu, muSun); % Sphere of influence considering the Sun as primary [km]
s.P = orbPeriod(s.sma, muSun); % Orbital period [s]
% Uranus
u.R = 25362.0; % Mean radius [km]
u.mu = 5.7940e6; % Gravitational parameter [kmˆ3/sˆ2]
u.sma = 2867.043e6; % Semi-major axis [km]
u.soi = soi(u.sma, u.mu, muSun); % Sphere of influence considering the Sun as primary [km]
u.P = orbPeriod(u.sma, muSun); % Orbital period [s]
% Neptune
n.R = 24622.0; % Mean radius [km]
n.mu = 6.8351e6; % Gravitational parameter [kmˆ3/sˆ2]
n.sma = 4514.953e6; % Semi-major axis [km]
n.soi = soi(n.sma, n.mu, muSun); % Sphere of influence considering the Sun as primary [km]
n.P = orbPeriod(n.sma, muSun); % Orbital period [s]
% % Pluto and charon (poor pluto is not a planet)
% p.R = 1188.0; % Mean radius [km]
% p.mu = 0.000870e6; % Gravitational parameter [kmˆ3/sˆ2]
% p.sma = 5869.656e6; % Semi-major axis [km]
% p.soi = soi(p.sma, p.mu, muSun); % Sphere of influence considering the Sun as primary [km]
% p.P = orbPeriod(p.sma, muSun); % Orbital period [s]

% Planets
planets.me = me;
planets.v = v ;
planets.e = e ;
planets.ma = ma;
planets.j = j;
planets.s = s;
planets.u = u;
planets.n = n;
% planets.p = p;

% Assign one color to each planet
fn = fieldnames(planets);
map = jet(length(fn));
map = flip(map);
map(4,:) = map(3,:); % I want Mars to be dark orange
map(3,:) = [0 0 0]; % Earth is black
for ii=1:length(fn)
    planets.(fn{ii}).color = map(ii,:);
end

end

function output = soi(sma, mu, muPrim)
%SOI computes sphere of influence
output = sma*(mu/muPrim)^(2/5);

end

function chiPoly = getChiPoly()
%GETCHIPOLY Returns the coefficients of the 8th degree polynomial that best fits the mapping function
chiPoly.p1 = -1.900257328955304; % [1/rad]
chiPoly.p2 = 13.641336822986435; % [1/rad]
chiPoly.p3 = -41.163077357520052; % [1/rad]
chiPoly.p4 = 67.643845254554719; % [1/rad]
chiPoly.p5 = -65.428068498802816; % [1/rad]
chiPoly.p6 = 37.880633918965955; % [1/rad]
chiPoly.p7 = -12.450354159214706; % [1/rad]
chiPoly.p8 = 0.083872533726236; % [1/rad]
chiPoly.p9 = 2.215294268840057; % [1/rad]

end

function sails = getSailsParams(srp)
%GETSAILSPARAMS Returns the a_0 and sigma of the flown solar sail missions to date
% IKAROS
ikaros.load = 310/(196*1e-6); % [kg/km^2]
ikaros.a0 =  load2a0(ikaros.load, srp); % [km/s^2]
% NanoSail-D2 (source: https://en.wikipedia.org/wiki/NanoSail-D2)
nanosaild.load = 4/(10*1e-6); % [kg/km^2]
nanosaild.a0 =  load2a0(nanosaild.load, srp); % [km/s^2]
% LightSail-2 (source: https://www.eoportal.org/satellite-missions/lightsail-2#spacecraft)
lightsail2.load = 5/(32*1e-6); % [kg/km^2]
lightsail2.a0 = load2a0(lightsail2.load, srp); % [km/s^2]
% LightSail-1 (assumed equal to LightSail-2)
lightsail1.a0 = lightsail2.a0; % [km/s^2]
lightsail1.load = lightsail2.load; % [kg/km^2]
% NEA scout (source: https://en.wikipedia.org/wiki/Near-Earth_Asteroid_Scout)
neascout.load = 14/(85*1e-6); % [kg/km^2]
neascout.a0 = load2a0(neascout.load, srp); % [km/s^2]
% ACS3 (source: data from NASA)
acs3.a0 = 0.045*1e-6; % [km/s^2]
acs3.load = a02load(acs3.a0, srp); % [kg/km^2]

% Output
sails.ikaros = ikaros;
sails.nanosaild = nanosaild;
sails.lightsail1 = lightsail1;
sails.lightsail2 = lightsail2;
sails.neascout = neascout;
sails.acs3 = acs3;

end

function a0 = load2a0(sailLoad, srp)
%LOAD2A0 Transforms sail load in [kg/km^2] to characteristic acceleration in [km/s^2]
a0 = 2*srp./sailLoad;

end

function sailLoad = a02load(a0, srp)
%A02LOAD Transforms characteristic acceleration in [km/s^2] to sail load in [kg/km^2]
sailLoad = 2*srp./a0;

end

function drhoPoly = getDrhoPoly()
%GETDRHOPOLY Returns the coefficients of the various polynomials that constitute the eta(dnu) envelope
drhoPoly.etaUpBest = [-0.061225179537496   0.281407601566408   0.366731636905808   0.000796829443216];
drhoPoly.etaUpWorst = [-6.810869349822958  17.886026812362250 -17.366950927475880   7.461553959036435  -0.736603745237182   0.159350171039189   0.013574703976377   0.000007739985217];
drhoPoly.etaDownBest = [0.019140172933931   0.062149393183795  -0.395940075916376   0.000403101999396];
drhoPoly.etaDownWorst = [0.750038253960143  -2.001360303317800   2.505544660451870  -1.611522552164965   0.054811325393119  -0.015979686340238  -0.000094807844629];
drhoPoly.G0 = 0.513195831412861; % ACS3 GEO (r0 = 42164 [km], a0 = 0.05 [mm/s^2])

end