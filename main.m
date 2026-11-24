% -------------- Enhanced Solar-sail Q-law  --------------------%
%
% --------------------------------------------------------------%
clc
clear all
close all
format long

% ------------------------------------------------------------- %
%                Default properties                             %
% ------------------------------------------------------------- %
currentpath = pwd; addpath([currentpath, '/..']);
addpath([currentpath, '/src/auxtools']);
addpath([currentpath, '/src/astrodynamics']);
addpath([currentpath, '/src/plotting']);
addpath([currentpath, '/src/solarsail']);  

% ------------------------------------------------------------- %
%                        Constants                              %
% ------------------------------------------------------------- %
const = getConstants();
ssprm = getSailParams();

% ------------------------------------------------------------- %
%                  Code controls                                %
% ------------------------------------------------------------- %
SimCase = 'testLeoSSO'; % Transfer in dawn-dusk SSO from 800 to 900 km
% SimCase = 'testLeoSSOfast'; % Transfer in dawn-dusk SSO from 800 to 810 km
SetTransferProb; 

plotResults = true;
saveOutput = false;

% ------------------------------------------------------------- %
%                Display header                                 %
% ------------------------------------------------------------- %
logicalStr = {'false', 'true'};
runTimeTag = datetime("now");
disp('% ------------------------------------------------------------- %');
disp('%                        main                                   %');
disp('% ------------------------------------------------------------- %');
disp(strcat("Executing from at time: ", char(runTimeTag,'yyyyMMdd-HHmmssS'))); 
disp(strcat("  ",pwd));
disp(" ");
disp("With code controls:");
disp(['  Simulation case: ', SimCase]);
disp(strcat("    -Save output: ",logicalStr{saveOutput + 1}));
disp("  Plotting:");
disp(strcat("    -Plot results: ",logicalStr{plotResults + 1}));
disp(" ");
disp('% ------------------------------------------------------------- %');

% ------------------------------------------------------------- %
%                  Dynamics settings                            %
% ------------------------------------------------------------- %
%% Transform initial COE to mEOE
mEOE0 = coe2mee(COE0)';
mEOEf = coe2mee(COEf)';

% ------------------------------------------------------------- %
%                     Propagate Target                          %
% ------------------------------------------------------------- %
%% Define gravity models
FM_PM    = FMmeoe_PointMass(const.mu); % Frameless class
FM_J2    = FMmeoe_OblateEarth_E(const.mu, const.Re, const.J2); % Computes the J2 perturbation in the Local frame if the mEOE are in Earth Frame (E)

if propTarget
    % Instantiate target's dynamics
    myForceModels = {FM_PM, FM_J2};
    targetDynamics = DynamicsOrbElem(const.mu, myForceModels, PropagationCoordinate, scaleFactor);
    
    % Propagate
    disp("Propagating target...")
    tic
    [~, ~, sol] = targetDynamics.PropagateSC_TL(Lspan(1) + 2*pi*[0 , 4], mEOEf, options, odeSolver);
    runtime = toc; % [s]
    fprintf('Elapsed time is %1.6f seconds \n\n', runtime);
    
    % Extract mean eccentricity from the solution
    N = length(sol.x);
    mEOEs = (sol.y(1:6,:).*[scaleFactor.l; ones(5,1)])';
    COEs = mee2coe(mEOEs);
    ecc = COEs(:,2);
    sma = COEs(:,1);
    argp = COEs(:,5);
    fprintf('mean ecc = %1.4e \n', mean(ecc))
    disp('% ------------------------------------------------------------- %');
    
    % Update final eccentricity
    COEf(2) = mean(ecc);
    mEOEf = coe2mee(COEf)';
    settings.COEf = COEf;

end

% ------------------------------------------------------------- %
%                     Define dynamics                           %
% ------------------------------------------------------------- %
%% Function handle that computes the rotation from Sun to Local frame
switch oeFrame
    case 'Ecliptic'
        f_QLS = @(t, COE) rotmatQLS(t, t0jd, COE); % For COEs in Ecliptic frame (P)
    case 'Earth'
        f_QLS = @(t, COE) rotmatQLS(t, t0jd, COE, 'Earth', const.obl); % For COEs in Earth frame (E)
    otherwise
        error('Frame not implemented')
end
QlawPrm.f_QLS = f_QLS;

% Function handle to function that computes the attitude of the sail that maximizes the total acceleration projection along the primer
switch SimCase
    case {'testLeoSSO','testLeoSSOfast'}
        f_maxAlongPrimer = @(t, coe, primer_L) optAttSrpOptOguri(t, coe, primer_L, QlawPrm.Cprms, f_QLS); % Optical SRP approximate solution, Oguri's approach
    otherwise
        error('Test case not implemented')
end

% Function handle to function that computes the primer vector direction:
f_cmptSCatt = @(t, mEOE, rkFlag) SRPctrl_LyaFB_mEOE(t, mEOE, mEOEf, f_maxAlongPrimer, QlawPrm, rkFlag); % Computes the attitude of the sail according to the Q-law. The input mEOE can be both in Earth (E) or Elciptic (P) frames

% Instantiate Solar sail model
FM_SRPECL     = FMmeoe_SRPctrl_Flat(QlawPrm, f_cmptSCatt, f_QLS, true); % Computes solar-sail acceleration in Local frame. The input mEOE can be both in Earth (E) or Elciptic (P) frames


% Selects force model depending on the simulation case. WARNING: Sail model always at the end for auxdata!
switch SimCase
    case {'testLeoSSO','testLeoSSOfast'}
        myForceModels = {FM_PM, FM_J2, FM_SRPECL};
    otherwise
        error('Test case not implemented')
end

% Selects convergence function depending on the simulation case
switch SimCase
    case {'testLeoSSO','testLeoSSOfast'}
        f_cmptConvergence = @(L, mEOE) cmptTargetErrorLeo(L, mEOE, COE0, COEf, const.mu, relconv);
    otherwise
        error('Convergence function not implemented')
end

% ------------------------------------------------------------- %
%               Display all relevant settings                   %
% ------------------------------------------------------------- %
%%
disp_settings;

% ------------------------------------------------------------- %
%                     Propagate Chaser                          %
% ------------------------------------------------------------- %
%% Instantiates the dynamics class
chaserDynamics = DynamicsOrbElem(const.mu, myForceModels, PropagationCoordinate, scaleFactor);

% Add convergence to target orbit event function
eventFunc = chaserDynamics.setEventFunc(struct('checkConvergence', 1, 'f_cmptConvergence', f_cmptConvergence, 'threshold', threshold));
options.Events = @(t,X) eventFunc(t,X);

% Integration in true longitude L = raan + argp + ta
disp("Propagating chaser...")
tic
[~, ~, sol] = chaserDynamics.PropagateSC_TL(Lspan, mEOE0, options, odeSolver);
% [~, ~, sol] = chaserDynamics.PropagateSC(timeSpan, mEOE0, options, odeSolver);
runtime = toc; % [s]
fprintf('Elapsed time is %1.6f seconds \n\n', runtime);
disp('% ------------------------------------------------------------- %');

% ------------------------------------------------------------- %
%           Evaluate accelerations over the solution            %
% ------------------------------------------------------------- %
%%
disp("Executing evalacc...")
auxdata.scaleFactor = scaleFactor;
auxdata.FM_SRP = myForceModels{end};
auxdata.QlawPrm = QlawPrm;
auxdata.Cprms = Cprms;
auxdata.a0 = a0;
auxdata.const = const;
auxdata.frame = oeFrame;
auxdata.t0jd = t0jd;
q = qsol2acc(sol, auxdata);
disp(" ")
disp('% ------------------------------------------------------------- %');

% ------------------------------------------------------------- %
% Store relevant information in the q structure before saving   %
% ------------------------------------------------------------- %
q.runtime = runtime;
q.drTarg = drTarg;
q.ssprm = ssprm;
q.t0jd = t0jd;
q.settings = settings;
q.settings.forceModels = myForceModels;

% ------------------------------------------------------------- %
%                        Save output                            %
% ------------------------------------------------------------- %
%% Execute if activated
if saveOutput
    % Create folder if it doesn't exist
    mkdir output

    % Save output
    h0 = sma0 - const.Re;
    fileName = ['qlaw_',num2str(a0*1e6),'a0_',num2str(h0),'h0_',num2str(drTarg),'dr_',num2str(rad2deg(COE0(3)),'%1.2f'),'inc_',num2str(rad2deg(COE0(4)),'%1.2f'),'raan_',char(runTimeTag,'yyyyMMdd-HHmmssS'),'.mat'];
    disp(['Saving output as ', fileName,'...']);
    save(fullfile('output', fileName));

    % End code execution if no plotting is required
    if ~plotResults
        return
    end

    disp(' ');
    disp('% ------------------------------------------------------------- %');

end

% ------------------------------------------------------------- %
%              Compute quantities for plotting                  %
% ------------------------------------------------------------- %
%%
disp("Executing qsol2plot...")
q = qsol2plot(q, auxdata);
disp(" ")
disp('% ------------------------------------------------------------- %');

% ------------------------------------------------------------- %
%                        CMD outputs                            %
% ------------------------------------------------------------- %
%% Time of flight
if q.converged
    disp(['ToF: ', num2str(max(abs(q.time))/const.sday), ' [day]'])
    disp(['revs: ', num2str(q.timerevs(end)), ' [revs]'])
    disp(" ")
else
    disp('Divergence!')
    disp(" ")
end
disp('% ------------------------------------------------------------- %');

% ------------------------------------------------------------- %
%                           Test                                %
% ------------------------------------------------------------- %
%%
tests

% ------------------------------------------------------------- %
%                           Plots                               %
% ------------------------------------------------------------- %
%%
if plotResults
    mainPlot
end