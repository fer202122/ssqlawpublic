function oenndot = safeguard(oenndot, COEtarg, prm)
%SAFEGUARD 
smaTarg = COEtarg(1);
smlVal = 1E-05 * prm.g0 * sqrt(smaTarg/prm.mu);
if any(oenndot > -smlVal)
    % disp('Applying safeguard')
    oenndot(oenndot > -smlVal) = -smlVal;
end

end