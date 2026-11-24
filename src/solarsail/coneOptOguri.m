function alphaOpt = coneOptOguri(phi, params)
%CONEOPTOPT Computes the normal vector cone angle that maximizes the 
% optical ss-acc along primer vector defined by the primer vector cone 
% angle phi. Problem solved following Oguri's approach

% Compute ideal optimum
alphaOptIdeal = coneOpt(phi);

% At zero the correction is not defined
if phi > 0 
    % Evaluate partial derivatives of J = a_s*lambda
    F = dJdcone(alphaOptIdeal, phi, params);
    F_alpha =  dJdcone2(alphaOptIdeal, phi, params);
    
    % Apply approximate correction to the ideal optimum
    alphaOpt = alphaOptIdeal - F./F_alpha;

else
    alphaOpt = alphaOptIdeal; 

end

end