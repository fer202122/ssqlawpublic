%% plotKeplerianOrbit: function description
function [fig, figprop, orbInfo, ah] = plotKeplerianOrbit(bodyName_or_COE, center, frame, et, fig, figprop, figsettings)
if ~exist('fig', 'var') || isempty(fig)
	fig = figure();
else
	figure(fig);
	hold on
end
if ~exist('figprop', 'var') || isempty(figprop)
	figprop = struct('Flag_3D', 1, 'tight', 1, 'skewAngle', [-50, 20]);
end
if ~exist('figsettings', 'var') || isempty(figsettings)
	figsettings.orbColor = 'b--';
else
	if ~isfield(figsettings, 'orbColor')
		figsettings.orbColor = 'b--';
	end
	if ~isfield(figsettings, 'LineWidth')
		figsettings.LineWidth = 3.0;
	end
end

if isnumeric(bodyName_or_COE)
	COE = bodyName_or_COE;
	sma = COE(1);
	ecc = COE(2);
	inc = COE(3);
	omg = COE(5);
	Omg = COE(4);
	TA0 = COE(6);
elseif isstruct(bodyName_or_COE)
	COE = bodyName_or_COE;
	sma = COE.a;
	ecc = COE.e;
	inc = COE.i;
	omg = COE.o;
	Omg = COE.O;
	TA0 = COE.TA;
elseif ischar(bodyName_or_COE)
	if ~exist('frame', 'var') || isempty(frame)
		frame = 'ECLIPJ2000';
	end
	if ~exist('et', 'var') || isempty(et)
		et0 = cspice_str2et('Jan 1 2000');
		xCartesian = cspice_spkezr(bodyName_or_COE, et0, frame, 'none', center);
	else
		xCartesian = cspice_spkezr(bodyName_or_COE, et, frame, 'none', center);
		orbit_rvec = xCartesian(1:3,:);
		orbit_rs = norm_vecs(orbit_rvec,1);
	end
	GM = cspice_bodvrd(center, 'GM', 1);
	[sma,ecc,~,TA0,inc,omg,Omg,~,~,~,~] = rv2elem(xCartesian(1:3,:),xCartesian(4:6,:),GM);
end

if ~exist('orbit_rvec', 'var')
	if ~exist('et', 'var') || isempty(et) || numel(et)==1
		et = nan;
		orbitTA = linspace(0,2*pi, 200);
	else
		GM = cspice_bodvrd(center, 'GM', 1);
		n_orbit = sqrt(GM/sma^3);
		orbitTA = TA0 + (et - et(1))*n_orbit;
	end
	orbit_rs = sma*(1-ecc^2)./(1+ecc*cos(orbitTA));
	RotOrb = Rot_space(3, Omg) * Rot_space(1, inc) * Rot_space(3,omg);
	ehat = RotOrb*[1;0;0];
	eper = RotOrb*[0;1;0];
	orbit_rvec = orbit_rs.*(cos(orbitTA).*ehat + sin(orbitTA).*eper);
end

ah = plot3(orbit_rvec(1,:),orbit_rvec(2,:),orbit_rvec(3,:), figsettings.orbColor, 'LineWidth', figsettings.LineWidth);
hold off


orbInfo = struct('sma', sma, 'ecc', ecc, 'inc', inc, 'TA0', TA0, 'omg', omg, 'Omg', Omg, ...
				 'et', et, 'orbit_rvec', orbit_rvec, 'orbit_rs', orbit_rs);
end