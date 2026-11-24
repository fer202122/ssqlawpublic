%% Rot_body: rotation about body axis
function [rotmat] = Rot_body(axis, angle_rad)
if axis == 'x' || axis == 1
	rotmat = rotx(angle_rad);
elseif axis == 'y' || axis == 2
	rotmat = roty(angle_rad);
elseif axis == 'z' || axis == 3
	rotmat = rotz(angle_rad);
end
end

function rotmat = rotx(angle_rad)

% rotate in the direction of y->z, counter-clockwise
if numel(angle_rad) == 1
	rotmat = [	1, 0, 0; ...
				0, cos(angle_rad), sin(angle_rad); ...
				0, -sin(angle_rad), cos(angle_rad)];
else
	N = numel(angle_rad);
	angle_rad = angle_rad(:)';
	rotmat_temp = [ones(1,N); zeros(1,N); zeros(1,N); zeros(1,N); cos(angle_rad); -sin(angle_rad); zeros(1,N); sin(angle_rad); cos(angle_rad)];
	rotmat = reshape(rotmat_temp, 3,3,N);
end
end

function rotmat = roty(angle_rad)

% rotate in the direction of z->x, counter-clockwise
if numel(angle_rad) == 1
	rotmat = [	cos(angle_rad), 0, -sin(angle_rad); ...
				0, 1, 0; ...
				sin(angle_rad), 0, cos(angle_rad)];
else
	N = numel(angle_rad);
	angle_rad = angle_rad(:)';
	rotmat_temp = [cos(angle_rad); zeros(1,N); sin(angle_rad); zeros(1,N); ones(1,N); zeros(1,N); -sin(angle_rad); zeros(1,N); cos(angle_rad)];
	rotmat = reshape(rotmat_temp, 3,3,N);
end
end

function rotmat = rotz(angle_rad)

% rotate in the direction of x->y, counter-clockwise
if numel(angle_rad) == 1
	rotmat = [	cos(angle_rad), sin(angle_rad), 0; ...
				-sin(angle_rad), cos(angle_rad), 0; ...
				0, 0, 1];
else
	N = numel(angle_rad);
	angle_rad = angle_rad(:)';
	rotmat_temp = [cos(angle_rad); -sin(angle_rad); zeros(1,N); sin(angle_rad); cos(angle_rad); zeros(1,N); zeros(1,N); zeros(1,N); ones(1,N)];
	rotmat = reshape(rotmat_temp, 3,3,N);
end
end