% ------------------------------------------------------------- %
%                    Simulation cases                           %
% ------------------------------------------------------------- %
switch SimCase
    % ------------------------------------------------------------- %
    %               Long transfer in dawn-dusk SSO                  %
    % ------------------------------------------------------------- %
    case 'testLeoSSO'
    PropagationCoordinate = 'mEOE'; 
    oeFrame = 'Earth';
    meanf = true;
    propTarget = true;

    % Sailcraft parameters
    Cprms = [1.476000000000000, 0.067932380952381, 0.262000000000000]; % [-] ACS3
    alpMinMax = [0, pi/2]; % [rad] Min/Max cone angle
    % ---- Characteristic acceleration ------------------------------------
    a0 = 0.045e-6; % Characteristic acceleration (ACS3) [km/s^2]

    % Initial orbit:
    % ---- SMA ------------------------------------------------------------
    sma0 = 800.0 + const.Re; % [km] Input orbital altitude
    % ---- ecc ------------------------------------------------------------
    ecc0 = 0.0; % [-]
    % ---- inc ------------------------------------------------------------
    coe.sma = sma0; coe.ecc = ecc0; inc0 = sunSincInc(coe); % [rad] Sun-syncrhonous inclination in Earth Frame (E)
    % ---- argp -----------------------------------------------------------
    argp0 = deg2rad(0); % [rad]
    % ---- raan -----------------------------------------------------------
    raan0 = deg2rad(90);
    % ---- tra ------------------------------------------------------------
    tra0 = deg2rad(45); % [rad] With 45 deg the initial SMA matches the average SMA
    % ---------------------------------------------------------------------
    COE0 = [sma0, ecc0, inc0, raan0, argp0, tra0]; % a, e, i, Omega, omega, TrueAnomaly

    % Target orbit:
    COEf = COE0;
    % ---- SMA ------------------------------------------------------------
    drTarg = 50; % [km]
    smaf = sma0 + drTarg; % [km] Add a drTarg to the initial SMA
    % ---- ecc ------------------------------------------------------------
    eccf = 0; % [-]
    % ---------------------------------------------------------------------
    COEf(1) = smaf; COEf(2) = eccf; % a, e, i, argp, raan

    % ---- initial epoch in julian date -----------------------------------
    thetas0 = 0; % [rad]
    t0jd = const.autEquiJd + thetas0/(2*pi)*const.syear/const.sday; % [jd]

	% Q-law settings:
    % ---- Weights per element --------------------------------------------
    Weights = [1; 1];
    % Q-law parameters
    const = getConstants();
	QlawPrm = struct(...
        'W', Weights,...                                       % Weights
        'g0', a0/2, 'Cprms', Cprms, 'alpMinMax', alpMinMax,... % Sail parameters
		'mu', const.mu,...                                     % Earth's gravitational parameter
        'COE0', COE0,...                                       % Initial COEs
        'const', const,...                                     % Physical constants
        't0jd', t0jd,...                                       % Initial time
        'meanf', meanf);                                       % Use rolling mean filter?                                          
    % Convergence criteria threshold (i.e., relative error)
    threshold = 1e-02; relconv = [1, 1]; % [-] Works always

    % Iintegration time-span:
    % Integration in true longitude ---------------------------------------
    revs = 700; % [revs] Number of orbital revolutions (LEO)
    T0 = 2*pi*sqrt(sma0^3/const.mu);
    t_end = T0*revs; % [s] Final time for scaling
    trueLon0 = wrapTo2Pi(sum(COE0(4:6)));
    Lspan = trueLon0 + 2*pi*[0 ,revs];

    % Integration settings:
    % Solver --------------------------------------------------------------
    odeSolver = @ode45;
    % Fixed step ----------------------------------------------------------
    Lstep = deg2rad(5)/(2*pi*revs);
    options = odeset('RelTol',1e+8,'AbsTol',1e+8, 'Events',@(xi, yi) DoNothing, 'OutputFcn', @odeplot, 'OutputSel', [1 2 3 7], 'MinStep', Lstep,'MaxStep', Lstep, 'InitialStep', Lstep);
    % Scaling -------------------------------------------------------------
    scaleFactor = struct('t', t_end, 'l', sma0,'L',2*pi*revs); % Like GPOPS scaling

    % ------------------------------------------------------------- %
    %               Short transfer in dawn-dusk SSO                 %
    % ------------------------------------------------------------- %
    case 'testLeoSSOfast'
    PropagationCoordinate = 'mEOE'; 
    oeFrame = 'Earth';
    meanf = true;
    propTarget = true;

    % Sailcraft parameters
    Cprms = [1.476000000000000, 0.067932380952381, 0.262000000000000]; % [-] ACS3
    alpMinMax = [0, pi/2]; % [rad] Min/Max cone angle
    % ---- Characteristic acceleration ------------------------------------
    a0 = 0.045e-6; % Characteristic acceleration (ACS3) [km/s^2]

    % Initial orbit:
    % ---- SMA ------------------------------------------------------------
    sma0 = 800.0 + const.Re; % [km] Input orbital altitude
    % ---- ecc ------------------------------------------------------------
    ecc0 = 0.0; % [-]
    % ---- inc ------------------------------------------------------------
    coe.sma = sma0; coe.ecc = ecc0; inc0 = sunSincInc(coe); % [rad] Sun-syncrhonous inclination in Earth Frame (E)
    % ---- argp -----------------------------------------------------------
    argp0 = deg2rad(0); % [rad]
    % ---- raan -----------------------------------------------------------
    raan0 = deg2rad(90);
    % ---- tra ------------------------------------------------------------
    tra0 = deg2rad(45); % [rad] With 45 deg the initial SMA matches the average SMA
    % ---------------------------------------------------------------------
    COE0 = [sma0, ecc0, inc0, raan0, argp0, tra0]; % a, e, i, Omega, omega, TrueAnomaly

    % Target orbit:
    COEf = COE0;
    % ---- SMA ------------------------------------------------------------
    drTarg = 10; % [km]
    smaf = sma0 + drTarg; % [km] Add a drTarg to the initial SMA
    % ---- ecc ------------------------------------------------------------
    eccf = 0; % [-]
    % ---------------------------------------------------------------------
    COEf(1) = smaf; COEf(2) = eccf; % a, e, i, argp, raan

    % ---- initial epoch in julian date -----------------------------------
    thetas0 = 0; % [rad]
    t0jd = const.autEquiJd + thetas0/(2*pi)*const.syear/const.sday; % [jd]

	% Q-law settings:
    % ---- Weights per element --------------------------------------------
    Weights = [1; 1];
    % Q-law parameters
    const = getConstants();
	QlawPrm = struct(...
        'W', Weights,...                                       % Weights
        'F', a0*cos(const.coneOpt)^2,...                       % Ion engine equivalent
        'g0', a0/2, 'Cprms', Cprms, 'alpMinMax', alpMinMax,... % Sail parameters
		'mu', const.mu,...                                     % Earth's gravitational parameter
        'COE0', COE0,...                                       % Initial COEs
        'const', const,...                                     % Physical constants
        't0jd', t0jd,...                                       % Initial time
        'meanf', meanf);                                       % Use rolling mean filter?                                          
    % Convergence criteria threshold (i.e., relative error)
    threshold = 1e-02; relconv = [1, 1]; % [-] Works always

    % Iintegration time-span:
    % Integration in true longitude ---------------------------------------
    revs = 700; % [revs] Number of orbital revolutions (LEO)
    T0 = 2*pi*sqrt(sma0^3/const.mu);
    t_end = T0*revs; % [s] Final time for scaling
    trueLon0 = wrapTo2Pi(sum(COE0(4:6)));
    Lspan = trueLon0 + 2*pi*[0 ,revs];

    % Integration settings:
    % Solver --------------------------------------------------------------
    odeSolver = @ode45;
    % Fixed step ----------------------------------------------------------
    Lstep = deg2rad(5)/(2*pi*revs);
    options = odeset('RelTol',1e+8,'AbsTol',1e+8, 'Events',@(xi, yi) DoNothing, 'OutputFcn', @odeplot, 'OutputSel', [1 2 3 7], 'MinStep', Lstep,'MaxStep', Lstep, 'InitialStep', Lstep);
    % Scaling -------------------------------------------------------------
    scaleFactor = struct('t', t_end, 'l', sma0,'L',2*pi*revs); % Like GPOPS scaling
    
    otherwise
        error('Case not implemented')

end


% ------------------------------------------------------------- %
%          Store all settings in a structure                    %
% ------------------------------------------------------------- %
%% Get a list of all variables
vars = whos;            % 
settings = struct;     % Initialize an empty structure

for i = 1:length(vars)
    varName = vars(i).name;
    settings.(varName) = eval(varName); % Use eval to access the variable's value
end
