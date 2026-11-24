function [u_C, angles] = bubbleContourOptical(aa, params)
%BUBBLECONTOUR Compute 3D control set
% Primer vector in Cranked Frame (C) (i.e., innerq grows counterclockwise from -y-axis)
N = 500;
innerq = linspace(0,pi,N)';
q_C = [sin(innerq), cos(innerq), zeros(N,1)];

% Primer vector in Inertial Frame (I)
QIC = rotmatQIC(aa);
q_I = matDotVec(QIC, q_C);

% Primer vector cone (psi) and clock angles in I. Note that primer clock angle is the same as the sail's clock angle. 
psi = acos(q_I(:,3));
clock = atan2(q_I(:,2), q_I(:,1));

% Sail's cone angle that yields maximum projection along primer in I
% cone = coneOpt(psi);
cone = coneOptOp(psi, params);

% Replace sail off with NaNs
cone(cone == pi/2) = NaN;

% as/a0 = u vector in I
% Optical parameters
C1 = params.C1;
C2 = params.C2;
C3 = params.C3;
C  = C1 + C2 + C3;
% Normal vector in I
nx_I = sin(cone).*cos(clock);
ny_I = sin(cone).*sin(clock);
nz_I = cos(cone);
% Control vector
u_I = zeros(N,3);
u_I(:,1) =  1/C*nz_I.*((C1*nz_I + C2).*nx_I     );
u_I(:,2) =  1/C*nz_I.*((C1*nz_I + C2).*ny_I     );
u_I(:,3) =  1/C*nz_I.*((C1*nz_I + C2).*nz_I + C3);

% n vector in C. Represents the control set that maximizes the Hamiltonian
QCI = pagetranspose(QIC);
u_C = matDotVec(QCI, u_I);
u_C = [u_C; [-1*u_C(end-1:-1:2,1), u_C(end-1:-1:2,2:3)]]; % Extend symmetrically excluding first and last

% Inner anlge and elevation of n
elevu = acos(u_C(:,3));           % Elevation angle gamma
inneru = atan2(u_C(:,1), u_C(:,2)); % Inner angle beta

% Interesting cone angles:
% Critical cone angle
const  = getConstants();
inneruAux = inneru(1:N);
coneC = spline(inneruAux(~isnan(inneruAux)), cone(~isnan(cone)), const.critcInner);
% Cone min/max
aa = wrapAa(aa);
coneMin = coneOptOp(pi/2 - aa, params);
coneMax = coneOptOp(pi/2 + aa, params);

% output
angles.innerq = [innerq; -1*innerq(end-1:-1:2)];
angles.inneru = inneru;
angles.psi = [psi; psi(end-1:-1:2)];
angles.elevu = elevu;
angles.coneMax = coneMax;
angles.cone = [cone; cone(end-1:-1:2)];
angles.coneMin = coneMin;
angles.coneC = coneC;

end