function [att, auxout] = no_qlaw
% no_qlaw sets all qlaw-related outputs to zero
att = zeros(2,1);
auxout.primer_L = zeros(3,1);
auxout.Q = 0.0;
auxout.dQdx = 0.0;
auxout.Bmat_L = zeros(2,2);
auxout.K = zeros(2,2);
auxout.tau = zeros(2,1);
auxout.OrbDist = zeros(2,1);
auxout.oenndot = zeros(2,1);
auxout.Lopt = zeros(2,1);
auxout.dtaudx = zeros(2,2);
auxout.oenndotAn = zeros(2,1);

end