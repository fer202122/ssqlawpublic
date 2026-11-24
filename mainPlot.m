% ------------ Plot Enhanced Solar-sail Q-law ------------------%
% 
% --------------------------------------------------------------%
close all

% ------------------------------------------------------------- %
%                Default properties                             %
% ------------------------------------------------------------- %
start_up; % Sets properties for the figures

% ------------------------------------------------------------- %
%                Display header                                 %
% ------------------------------------------------------------- %
disp('% ------------------------------------------------------------- %');
disp('%                     mainPlot                                  %');
disp('% ------------------------------------------------------------- %');

% ------------------------------------------------------------- %
%            Complete q structure if required                   %
% ------------------------------------------------------------- %
if ~isfield(q,'timerevs')
    disp("Executing qsol2plot...")
    q = qsol2plot(q,auxdata);
end

% ------------------------------------------------------------- %
%          Substitude zeros by NaNs during eclipse              %
% ------------------------------------------------------------- %
fns = fieldnames(q);
for ii=1:numel(fns)
    field = q.(fns{ii});
    sz = size(field);
    if sz(1) == q.N && ~strcmp(fns{ii},'eclipse') && ~strcmp(fns{ii},'contRevs_C') && ~strcmp(fns{ii},'contAngRevs') && ~strcmp(fns{ii},'r_S') && ~strcmp(fns{ii},'r_P') && ~strcmp(fns{ii},'r_E') && ~strcmp(fns{ii},'timerevs') && ~strcmp(fns{ii},'timeyears') && ~strcmp(fns{ii},'time') && ~strcmp(fns{ii},'timedays') && ~strcmp(fns{ii},'truelong') && ~strcmp(fns{ii},'COEs_E') && ~strcmp(fns{ii},'COEsm_E') && ~strcmp(fns{ii},'powSs')
        field(q.eclipse,:) = NaN;
        q.(fns{ii}) = field;

    end

end

%% Trajectory in Sun frame (S)
primerHat_S = makeUnit(q.primer_S); % Make primer unit vector

figs = {};
figprops = {};
elemColor = {'c', 'r', 'g'};
figs{numel(figs)+1} = figure();
figprops{length(figprops)+1} = struct('tight', 0, 'Flag_3D', 1, 'skewAngle', [-53, 13]);
    map = copper(q.N);
    hold on
    traj = plot3(q.r_S(:,1), q.r_S(:,2), q.r_S(:,3), 'k');
    ssa = quiver3(q.r_S(:,1),q.r_S(:,2),q.r_S(:,3),q.as_S(:,1),q.as_S(:,2),q.as_S(:,3),'DisplayName','Solar-sail acceleration','AutoScale','on','Color',map(end,:));
    prim = quiver3(q.r_S(:,1),q.r_S(:,2),q.r_S(:,3),primerHat_S(:,1),primerHat_S(:,2),primerHat_S(:,3),'DisplayName','Primer vector','AutoScale','on','Color',map(1,:));
    h0 = scatter3(q.r_S(1,1), q.r_S(1,2), q.r_S(1,3), 'b^', 'filled');
    hf = scatter3(q.r_S(end,1), q.r_S(end,2), q.r_S(end,3), 'r^', 'filled');
    earth_sphere(); alpha 0.5; hold on;
    scatter3(0,0,0, [], 'k', 'Filled')
    [~,~,~,p_orb0] = plotKeplerianOrbit(q.COEs_S(1,:), [], [], [], figs{numel(figs)}, figprops{numel(figprops)}, struct('orbColor', 'r', 'LineWidth', 1.5));
    [~,~,~,p_orbf] = plotKeplerianOrbit(q.COEs_S(end,:), [], [], [], figs{numel(figs)}, figprops{numel(figprops)}, struct('orbColor', 'g', 'LineWidth', 1.5));
    legend([p_orb0, traj, p_orbf, h0, hf, ssa, prim], {'Departure orbit', 'Transfer','Arrival orbit', 'Initial state', 'Final state','$a_s$','$p$'}, 'location', 'best')
    grid on
    xlabel('$x_s$, [km]')
    ylabel('$y_s$, [km]')
    zlabel('$z_s$, [km]')
    view([120, 20])
    axis equal
    axis tight

%% Evolution of the primer unit vector components in Local frame (L)
primerHat_L = makeUnit(q.primer_L);

figure()
    labels = ["$\hat{p}_{x_l}$ [-]","$\hat{p}_{y_l}$ [-]","$\hat{p}_{z_l}$ [-]"];
    for jj=1:length(labels)
        subplot(3, 1, jj);
            hold on
            h2 = plot(q.timeyears,  primerHat_L(:,jj),'Color','k');
            ylabel(labels(jj))
            hold off
    end
    arrayfun(@(x) axis(x,'tight'), findobj(gcf,'Type','axes'));
    arrayfun(@(x) grid(x,'on'), findobj(gcf,'Type','axes'));
    sgtitle('Primer components in Local frame');
    xlabel('Time [years]')

%% Evolution of semi-major axis, eccentricity, and aspect angle
figure();
    idxs = [1 2 99];
    element = ["$\Delta a$ [km]" "$e$ [-]" "$\beta$ [deg]"];
    for jj=1:length(idxs)
        subplot(length(idxs), 1, jj);
            hold on;
            idx = idxs(jj);
            if idx==1 % Delta SMA
                h1=yline(q.drTarg,'r--','LineWidth',2,'DisplayName','Target');
                plot(q.timeyears, q.COEs_P(:,idx) - q.COEs_P(1,idx),'Color','k','LineWidth',2)
                legend(h1,'Location','SouthEast')
            elseif idx==2 % Ecc
                plot(q.timeyears, q.COEs_P(:,idx),'Color','k','LineWidth',2)
                h2=yline(q.settings.COEf(2),'r--','LineWidth',2,'DisplayName','Target');
                legend(h2,'Location','SouthEast')
            elseif idx==99 % Aspect angle
                plot(q.timeyears, rad2deg(wrapAa(q.aa)),'Color','k','LineWidth',2)
            end
            ylabel(element(jj));
            hold off;
    end
    xlabel('Time [years]');
    arrayfun(@(x) axis(x,'tight'), findobj(gcf,'Type','axes'));
    arrayfun(@(x) grid(x,'on'), findobj(gcf,'Type','axes'));
    sgtitle('$a$, $e$, and $\beta$');

 %% Evolution of semi-major axis and eccentricity relative error
 sma = figure();
    idxs = [1 2];
    element = ["$\frac{\vert \bar{a} - \bar{a_t} \vert}{\vert \bar{a_0} - \bar{a_t} \vert}$"...
               "$\frac{\vert \bar{e} - \bar{e_t} \vert}{\vert \bar{e_0} - \bar{e_t} \vert}$"];
    for jj=1:length(idxs)
        subplot(length(idxs), 1, jj);
        hold on;
        idx = idxs(jj);
        if idx==1 % SMA relative error w.r.t Delta SMA
            plot(q.timeyears, abs((q.COEsm_E(:,idx) - q.settings.COEf(idx))/q.drTarg),'Color','k','LineWidth',1.5)
        elseif idx==2 % Ecc relative error
            plot(q.timeyears, abs((q.COEsm_E(:,idx) - q.settings.COEf(idx))/q.settings.COEf(idx)),'Color','k','LineWidth',1.5)
        end
        hold on
        h = yline(q.settings.threshold*q.settings.relconv(idx),'r--','LineWidth',1.5,'DisplayName','Threshold');
        hold off
        yscale log
        ylabel(element(jj),'Interpreter','Latex');
        hold off;
        legend(h,'Location','Best')
    end
    xlabel('Time [years]');
    arrayfun(@(x) axis(x,'tight'), findobj(gcf,'Type','axes'));
    arrayfun(@(x) grid(x,'on'), findobj(gcf,'Type','axes'));
    sgtitle(['$a$ and $e$ relative error, $e_0$=0 [-], ','$e_t$=', num2str(q.settings.COEf(2),'%1.3e'),' [-]']);

%% Classical orbital elements in Earth frame
figure()
    element = ["$a$ [km]" "$e$ [-]" "$i$ [deg]" "$\Omega$ [deg]" "$\omega$ [deg]"];
    COEdeg = q.COEs_E;
    COEdeg = [COEdeg(:,1:2), rad2deg(COEdeg(:,3:5))];
    COEmdeg = q.COEsm_E;
    COEmdeg = [COEmdeg(:,1:2), rad2deg(COEmdeg(:,3:5))];
    for jj=1:length(element)
        subplot(length(element), 1, jj);
            h1=plot(q.timeyears, COEdeg(:,jj),'Color','k','DisplayName','Osculating');
            hold on
            h2=plot(q.timeyears, COEmdeg(:,jj),'Color','b','LineWidth',1.5,'DisplayName','Average');
            hold off
            ylabel(element(jj));
            if jj==1
                legend([h1,h2],'Location','SouthEast')
            end
    end
    xlabel('Time [years]');
    arrayfun(@(x) axis(x,'tight'), findobj(gcf,'Type','axes'));
    arrayfun(@(x) grid(x,'on'), findobj(gcf,'Type','axes'));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.45, 0.9]*0.9);
    sgtitle('Classical orbital elements in Earth frame');

%% Control cone and clock angles in Sun frame
plotControlHistory(q.timeyears, (q.att'), struct('controlType', 'solarSail', 't_unit', '[years]'));
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.45, 0.5]*0.9);
xlabel('Time [years]')
