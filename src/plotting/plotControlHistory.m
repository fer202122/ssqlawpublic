function [figs, figprops] = plotControlHistory(t_his, u_his, option)
if ~isfield(option, 'controlType')
	option.controlType = 'lowThrust';
end

switch option.controlType
	case 'lowThrust'
		if ~isfield(option, 'Nu')
			Nu = 3;
		end
		if ~isfield(option, 'labelList')
			option.labelList = {	'$$a_1$$', ...
									'$$a_2$$', ...
									'$$a_3$$'};
		end
	case 'solarSail'
		if ~isfield(option, 'Nu')
			Nu = 2;
		end
		option.labelList = {	'$$\alpha\ \mathrm{[deg]}$$', ...
								'$$\beta\ \mathrm{[deg]}$$'};
		u_his(1,:) = rad2deg(u_his(1,:));
		u_his(2,:) = rad2deg(wrapToPi(u_his(2,:)));
end


figs{1} = figure();
figprops{1} = struct('tight', 0, 'Flag_3D', 0);
hold on
for ii=1:Nu
	subplot(Nu,1,ii);
	hold on
	pnom = plot(t_his, u_his(ii,:), 'k-');
	ylabel(option.labelList{ii});
	grid on
	xlim([t_his(1), t_his(end)])
	if ii==Nu
		xlabel(['$$t$$, ', option.t_unit]);
	end
end
axis tight



end