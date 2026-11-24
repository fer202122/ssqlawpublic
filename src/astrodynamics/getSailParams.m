function output = getSailParams()
%GETSAILPARAMS Group sail parameters under the same function
    % ACS3 moments of inertia
    Ixx =  21.0; % [kg*m^2]
    Iyy =  10.5; % [kg*m^2]
    Izz =  10.5; % [kg*m^2]
    % Drag model related parameters
    sigT = 0.8; % Tangential momentum accomodation coefficient [-]
    sigN = 0.8; % Normal momentum accomodation coefficient [-]
    VR = 0.05;  % Ratio between the thermal motion of the molecules to the characteristic velocity of the spacecraft [-]
    eps = 1e-4; % Small constant to render CD and CLsimpl derivatives continuous (i.e., abs(x) ~ sqrt(eps + x^2))

    % Optical coefficients
    % % NEA scout
    % C1 = 1.711;
    % C2 = 0.002;
    % C3 = 0.145;
    % % Ideal
    % C1 = 2;
    % C2 = 0;
    % C3 = 0;
    % ACS3
    C1 = 1.476000000000000; % 2*0.9*0.82
    C2 = 0.067932380952381; % 0.79*0.9*(1-0.82) + (1 - 0.9)*(0.79*0.03 - 0.67*0.6)/(0.03+0.6)
    C3 = 0.262; %  1 - 0.9*0.82

    % Critical angles
    [incc, betac] = incCrit([C1, C2, C3]);
    [conec, phic] = coneCrit([C1, C2, C3]);
    output.incc = incc; 
    output.betac = betac;  
    output.conec = conec;  
    output.phic = phic; 

    % output
    output.Ixx = Ixx;
    output.Iyy = Iyy;
    output.Izz = Izz;
    output.sigT = sigT;
    output.sigN = sigN;
    output.VR = VR;
    output.eps = eps;
    output.C1 = C1;
    output.C2 = C2;
    output.C3 = C3;

end

