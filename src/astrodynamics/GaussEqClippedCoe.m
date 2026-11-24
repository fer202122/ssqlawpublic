function Bmat_L = GaussEqClippedCoe(GM, COE)
%GaussEqCoe Computes the non-singular Gauss dynamics matrix in COE

% Safeguard for eccentricity
if COE(2) < 1e-12
    % fprintf('Clipping eccentricity, sma = %1.2e \n', COE(1))
    COE(2) = 1e-12;
end

% Safeguard for inclination
if COE(3) < 1e-12
    % disp('Clipping inclination')
    COE(3) = 1e-12;
end

% Compute B matrix
Bmat_L = GaussPlanetaryEqs1(GM, COE);

end