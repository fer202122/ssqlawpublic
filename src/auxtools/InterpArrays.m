%% InterpArrays: function description
function [x_s] = InterpArrays(t, x, t_s, alg)
% x: time-profile of vector
if ~exist('alg', 'var')
	alg = 'interp1';
end
transposed = 0;
if numel(x(:,1))~=numel(t)
	transposed = 1;
	x = x';
	t = t';
	t_s = t_s';
end

%% Setting initial condition
tx = [t, x];
if t(1) < t(end)
	tx_ = [tx(find(diff(tx(:,1))>0),:); tx(end,:)];
else
	tx_ = [tx(1,:); tx(find(diff(tx(:,1))<0),:)];
end
x_s = zeros(numel(t_s),numel(x(1,:)));

switch alg
	case 'spline'
	for jj = 1:numel(x(1,:))
		x_s(:,jj) = ppval(spline(tx_(:,1),tx_(:,jj+1)),t_s);
	end

	case 'interp1'
	for jj = 1:numel(x(1,:))
		x_s(:,jj) = interp1(tx_(:,1),tx_(:,jj+1),t_s);
	end
end

if transposed
	x_s = x_s';
end

end
