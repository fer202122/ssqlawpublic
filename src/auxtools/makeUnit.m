function vecHat = makeUnit(vec)
    vecHat = vec./vecnorm(vec')';
end

