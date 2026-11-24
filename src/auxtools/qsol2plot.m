function q = qsol2plot(q, auxdata)
%% Unfold auxdata
QlawPrm = auxdata.QlawPrm;
Cprms = auxdata.Cprms;
mu = QlawPrm.mu;
const = auxdata.const;
frame = auxdata.frame;
t0jd = auxdata.t0jd;

% ----------------------------------------------------------------------- %
% Clasical orbital elements in various frames
% ----------------------------------------------------------------------- %
% Transform Modified Equinoctial Elements to Classical Orbital Elements
q.COEs_P = mee2coe(q.mEOEs_P);
q.COEs_E = mee2coe(q.mEOEs_E);

% Compute mean classical orbital elements
COEsm_E = zeros(size(q.COEs_E));
for ii=1:size(q.COEs_E,1)
    COEsm_E(ii,:) = rollingMean(q.time(ii), q.COEs_E(ii,:)', mu);
end
COEsm_P = zeros(size(q.COEs_P));
for ii=1:size(q.COEs_E,1)
    COEsm_P(ii,:) = rollingMean(q.time(ii), q.COEs_P(ii,:)', mu);
end
q.COEsm_E = COEsm_E;
q.COEsm_P = COEsm_P;

% ----------------------------------------------------------------------- %
% Position and velocity in various frames
% ----------------------------------------------------------------------- %
% Transform general orbital elements to position and velocity history
auxdata.mu = mu;
cart_P = mee2cart(q.mEOEs_P, auxdata);
q.r_P = cart_P(:,1:3);
q.v_P = cart_P(:,4:6);

% Position and velocity in Sun frame (S)
QPS = rotmatQPS(q.time, t0jd);
QSP = pagetranspose(QPS);
q.r_S = matDotVec(QSP, q.r_P);
q.v_S = matDotVec(QSP, q.v_P);

% Position and velocity in Earth frame (E)
QEP = rotmatQEP;
q.r_E = matDotVec(QEP, q.r_P);
q.v_E = matDotVec(QEP, q.v_P);

% ----------------------------------------------------------------------- %
% Motion of the Sun in various frames
% ----------------------------------------------------------------------- %
q.thetas = sunMotion(q.time, t0jd);
q.sunHat_P = [cos(q.thetas), sin(q.thetas), zeros(q.N,1)];

coe_P.inc = q.COEs_P(:,3);
coe_P.raan = q.COEs_P(:,4);
QCP = rotmatQCP(coe_P);
sunHat_C = matDotVec(QCP, q.sunHat_P);

q.sunHat_E = matDotVec(QEP, q.sunHat_P);

% Get COEs in Sun frame
q.COEs_S = q.COEs_P;
q.COEs_S(:,4) = q.COEs_S(:,4) - q.thetas;

% ----------------------------------------------------------------------- %
% Argument of eclipse in Ecliptic frame (P) / Cranked frame (C)
% ----------------------------------------------------------------------- %
% Obtain argument of latitude from the argument of perigee and the true anomaly
argl = q.COEs_P(:,5) + q.COEs_P(:,6);
argl = linearAngle(argl);

% Find the indices of the entry and exit points from eclipse
argEclipse = NaN(size(argl));
argEclipse(q.eclipse) = argl(q.eclipse);
firstArglIdxs = strfind(isnan(argEclipse)',[1 0])' + 1;
secondArglIdxs = strfind(isnan(argEclipse)',[0 1])';
if q.eclipse(1)
    firstArglIdxs = [1; firstArglIdxs];
end
if q.eclipse(end)
    secondArglIdxs = [secondArglIdxs; length(argl)];
end

% Compute "argument of eclipse" as the mean angular distance between entry and exit points
meanArgl = (argl(firstArglIdxs) + argl(secondArglIdxs))./2;

% Assign constan values (entry, mean, exit) to the entire eclipse window and wrap to 2*pi
argEclipse = [NaN(size(argl)), argEclipse, NaN(size(argl))];
for ii=1:length(meanArgl)
    argEclipse(firstArglIdxs(ii):secondArglIdxs(ii),1) = wrapTo2Pi(argl(firstArglIdxs(ii)));
    argEclipse(firstArglIdxs(ii):secondArglIdxs(ii),2) = wrapTo2Pi(meanArgl(ii));
    argEclipse(firstArglIdxs(ii):secondArglIdxs(ii),3) = wrapTo2Pi(argl(secondArglIdxs(ii)));
end
q.argEclipse = argEclipse;

% ----------------------------------------------------------------------- %
% Solar-sail acceleration in various frames
% ----------------------------------------------------------------------- %
% Rotate to Ecliptic (P) and Sun (S) frame 
switch frame
    case 'Earth'
        QLS = QlawPrm.f_QLS(q.time, q.COEs_E);
    case 'Ecliptic'
        QLS = QlawPrm.f_QLS(q.time, q.COEs_P);
    otherwise
        error('Frame not implemented');
end
QSL = pagetranspose(QLS);
QPL = pagemtimes(QPS,QSL);
as_S = matDotVec(QSL, q.as_L);
as_P = matDotVec(QPS, as_S);

% output
q.as_P = as_P;
q.as_S = as_S;
q.as_C = matDotVec(QCP, as_P);

% ----------------------------------------------------------------------- %
% Yaw and pitch angles of the SS-acceleration vector in Local frame
% ----------------------------------------------------------------------- %
q.yawss = atan2(q.as_L(:,1), q.as_L(:,2));
q.as_Luv = sqrt(q.as_L(:,1).^2 + q.as_L(:,2).^2);
q.pitchss = atan2(q.as_L(:,3), q.as_Luv);

% ----------------------------------------------------------------------- %
% Yaw and pitch angles of the angular momentum vector wrt the Sun vector
% ----------------------------------------------------------------------- %
% Angular momentum vector in Sun frame
hHat_S = vecCrossVec(q.r_S, q.v_S);
hHat_S = hHat_S./vecNorm(hHat_S);
hHatxy_S = vecNorm(hHat_S(:,[1,2]));

q.yawh = atan2(hHat_S(:,2), hHat_S(:,1));
q.pitchh = atan2(hHat_S(:,3), hHatxy_S);

% ----------------------------------------------------------------------- %
% Angular momentum vector in various frames
% ----------------------------------------------------------------------- %
q.hHat_S = hHat_S;
q.hHat_P = matDotVec(QPS, q.hHat_S);
q.hHat_E = matDotVec(QEP, q.hHat_P);

% ----------------------------------------------------------------------- %
% Times
% ----------------------------------------------------------------------- %
% Compute time in orbital revolutions
switch frame
    case 'Earth'
        if q.COEs_E(1,3) == 0
            argl = q.truelong;
            warning('We cannot use average J2 effect with zero inclination')
        else
            argl = q.truelong  - wrapTo2Pi(q.COEs_E(:,4)); 
        end

    case 'Ecliptic'
        if q.COEs_P(1,3) == 0
            argl = q.truelong;
        else
            argl = q.truelong  - wrapTo2Pi(q.COEs_P(:,4));
        end

    otherwise
        error('Frame not implemented')
end
% We assume the number of orbital revolutions is independent of the frame
q.timerevs = (argl - argl(1))/(2*pi);

% Compute time in days
q.timedays = q.time/const.sday;

% Time in years
q.timeyears = q.time/const.syear;

% ----------------------------------------------------------------------- %
% Primer vector in various frames
% ----------------------------------------------------------------------- %
q.primer_P = matDotVec(QPL, q.primer_L);
q.primer_S = matDotVec(QSP, q.primer_P);

% ----------------------------------------------------------------------- %
% Primer clock angle
% ----------------------------------------------------------------------- %
% Compute clock angle measured from z_s axis [-pi, pi]
q.pclock = atan2(q.primer_S(:,2), q.primer_S(:,3));

% ----------------------------------------------------------------------- %
% Velocity cone and clock angles
% ----------------------------------------------------------------------- %
vHat_S = makeUnit(q.v_S);
% Compute clock angle measured from z_s axis [-pi, pi]
q.vclock = atan2(vHat_S(:,2), vHat_S(:,3));
% Compute cone angle measured from z_s axis [0, pi/2]
q.vcone = acos(vHat_S(:,1));

% ----------------------------------------------------------------------- %
% Primer yaw angle
% ----------------------------------------------------------------------- %
q.yaw = atan2(q.primer_L(:,1), q.primer_L(:,2)); % [-pi, pi] Primer yaw angle measured from T axis

% ----------------------------------------------------------------------- %
% Aspect angle
% ----------------------------------------------------------------------- %
q.aa = aspectAngle([q.r_P, q.v_P], q.sunHat_P);

% ----------------------------------------------------------------------- %
% Projected contour of the force bubble in the orbital plane
% ----------------------------------------------------------------------- %
q.revs0 = floor(q.timerevs(end));
aaRevs = splitPerRev(q, 'aa');
q.sunHat_C = sunHat_C;
sunHatRevs_C = splitPerRev(q, 'sunHat_C');

params.C1 = Cprms(1);
params.C2 = Cprms(2);
params.C3 = Cprms(3);

contRevs_C = cell(q.revs0,1);
contAngRevs = cell(q.revs0,1);
for ii=1:q.revs0
    aaRev = aaRevs{ii};
    % [cont_B, contAng] = bubbleContour(aaRev(1));
    [cont_B, contAng] = bubbleContourOptical(aaRev(1), params);

    % Bubble angles
    contAngRevs{ii} = contAng;

    % Rotate to align with current Sunline
    sunHat_C = sunHatRevs_C{ii};
    angle = atan2(-sunHat_C(1,1), sunHat_C(1,2));
    QCB = [[  cos(angle) -sin(angle) 0];...
           [  sin(angle)  cos(angle) 0];...
           [ 0 0 1]];
    cont_C = matDotVec(QCB, cont_B);
    contRevs_C{ii} = cont_C;

end

% output
q.contRevs_C = contRevs_C;
q.contAngRevs = contAngRevs;

% ----------------------------------------------------------------------- %
% Osculating eccentricity vector and growth vector
% ----------------------------------------------------------------------- %
% Osculating eccentricity
h_P = vecCrossVec(q.r_P, q.v_P);
r = vecNorm(q.r_P);
eccVec_P = vecCrossVec(q.v_P, h_P)/mu - q.r_P./r;
q.eccVec_C = matDotVec(QCP,eccVec_P);

% Growth vector
q.growHat_C = -symAxis(q.sunHat_C);

% ----------------------------------------------------------------------- %
% SMA & eccentricity rates of change from primers
% ----------------------------------------------------------------------- %
apDot = zeros(size(q.time));
epDot = zeros(size(q.time));
for ii=1:length(q.time)
    COE_E = q.COEs_E(ii,:);

    B = GaussPlanetaryEqs1(const.mu, COE_E);

    primerSma_L = B(1,:); % Primer to change SMA
    primerEcc_L = B(2,:); % Primer to change eccenetricity

    apDot(ii) = vecDotVec(primerSma_L, q.as_L(ii,:));
    epDot(ii) = vecDotVec(primerEcc_L, q.as_L(ii,:));

end
q.apDot = apDot;
q.epDot = epDot;

end
