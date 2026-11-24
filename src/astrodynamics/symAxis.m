function symHat_C = symAxis(sunHat_C)
%SYMAXIS Computes the symmetry axis from the Sun vector in Cranked frame
N = length(sunHat_C(:,1));
symHat_C = vecCrossVec(sunHat_C, [zeros(N,1) zeros(N,1) ones(N,1)]);
symHatNorm = vecNorm(symHat_C);
if any(symHatNorm < 1e-6)
    warning('Cross product is almost zero!')
    symHatNorm(symHatNorm < 1e-6) = NaN;
end    
symHat_C = symHat_C./symHatNorm;

end