function q = qsol2acc(sol, auxdata)
%% Unfold auxdata
scaleFactor = auxdata.scaleFactor;
FM_SRP = auxdata.FM_SRP;
QlawPrm = auxdata.QlawPrm;
Cprms = auxdata.Cprms;
mu = QlawPrm.mu;
frame = auxdata.frame;
t0jd = auxdata.t0jd;

% output
q.a0 = auxdata.a0;
q.QlawPrm = QlawPrm;
q.Cprms = Cprms;
q.t0jd = t0jd;

%% Interpolate solution
% N = 40000;
% q.truelong = linspace( Lspan(1), Lspan(end), N )./scaleFactor.L;
% [q.mEOEs_P, q.time] = myDynamics.getTrajectory(sol, q.truelong);
% q.truelong = q.truelong*scaleFactor.L;
% N = length(q.time);

% Bypassing the getTrajectory function (because it does the wierd thing with the final time):
% Integration in True Longitude
N = length(sol.x);
q.N = N;
q.truelong = (sol.x*scaleFactor.L)';
q.time = (sol.y(7,:)*scaleFactor.t)';
mEOEs = (sol.y(1:6,:).*[scaleFactor.l; ones(5,1)])';
% % Integration in time
% N = length(sol.x);
% N = N;
% q.time = sol.x*scaleFactor.t;
% q.truelong = sol.y(6,:);
% q.mEOEs_P = sol.y(1:6,:).*[scaleFactor.l; ones(5,1)];

% ----------------------------------------------------------------------- %
% Convergence flag
% ----------------------------------------------------------------------- %
if isempty(sol.ie)
    q.converged = false;
else
    q.converged = true;
end

% ----------------------------------------------------------------------- %
% Modified Equinoctial Elements in various frames
% ----------------------------------------------------------------------- %
switch frame
    case 'Earth'
        % Transform mEOEs from Earth (E) to Ecliptic (P) Frame
        q.mEOEs_E = mEOEs;
        q.mEOEs_P = mEOEsE2P(mEOEs, mu);

    case 'Ecliptic'
        q.mEOEs_P = mEOEs;
        q.mEOEs_E = mEOEsP2E(mEOEs, mu);

    otherwise
        error('Frame not implemented')
end

% ----------------------------------------------------------------------- %
% Evaluate control dependent acceleration terms (ssacc and aeroacc)
% ----------------------------------------------------------------------- %
% Acceleration models related outputs
q.at_L = zeros(N,3); % Total acceleration
q.as_L = zeros(N,3); % Solar-sail acceleration in Local frame (L)
q.eclipse = false(N,1); % Eclipse logical variable
% Q-law related outputs
q.t_Q_Qdot = zeros(N,3); % Array of time, Q and Qdot
q.Ssma = zeros(N,1); % Value of the S function for the SMA
q.L_plot = zeros(N,1); % True longitude
q.dQdx = zeros(N,6); % Gradient of the Q function
q.primer_L = zeros(N,3); % Primer vector in Local frame (L)
q.att = zeros(N,2); % Attitude history vector
q.OEnndot = zeros(N,2); % Double minimization q per element
q.tau = zeros(N,2); % Double minimization q per element
q.Lopt = zeros(N,2); % Optimal true longitudes
q.OrbDist = zeros(N,2); % Optimal true longitudes
q.Bmat_L = zeros(2,2,N); % Optimal true longitudes
q.weights = zeros(N,2); % Used weights though matrix K
q.dtaudx = zeros(2,2,N); % tau partials
q.oenndotAn = zeros(N,2); % Double minimization q per element analytical
for ii=1:N
    % Compute ssa from the model used inside the propagation in Local frame (L)
    switch frame
    case 'Earth'
        mEOEi = q.mEOEs_E(ii,:)';
    case 'Ecliptic'
        mEOEi = q.mEOEs_P(ii,:)';
    otherwise
        error('Frame not implemented');
    end
  
    [q.at_L(ii,:), auxout] = FM_SRP.cmpt_accel(q.time(ii), mEOEi, false);
    % Acceleration models related outputs
    q.as_L(ii,:) = auxout.as_L;
    q.eclipse(ii) = auxout.eclipse;
    % Q-law related outputs
    q.att(ii,:) = auxout.att;
    q.primer_L(ii,:) = auxout.primer_L;
    q.t_Q_Qdot(ii,2) = auxout.Q;
    q.dQdx(ii,:) = auxout.dQdx;
    q.Ssma(ii,1) = 1;
    q.weights(ii,:) = [auxout.K(1,1), auxout.K(2,2)];
    q.tau(ii,:) = auxout.tau;
    q.OrbDist(ii,:) = auxout.OrbDist;
    q.OEnndot(ii,:) = auxout.oenndot;
    q.Lopt(ii,:) = auxout.Lopt;
    q.Bmat_L(:,:,ii) = auxout.Bmat_L(1:2,1:2);
    q.dtaudx(:,:,ii) = auxout.dtaudx;
    q.oenndotAn(ii,:) = auxout.oenndotAn;
    dQdt = -1*q.primer_L(ii,:)*q.as_L(ii,:)';
    q.t_Q_Qdot(ii,3) = dQdt;

end

end
