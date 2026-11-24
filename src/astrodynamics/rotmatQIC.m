function QIC = rotmatQIC(inc)
%ROTMATQIC Rotation matrix from Cranked to Inertial frame. This matrix is
% intended to be used for the first paper, when the Sun doesn't move and
% it's always in the -z direction. Not to be confused with the matrices 
% QPC etc

% Rotation around the x-axis
N = length(inc);
QIC = cat(3,[ ones(N,1), zeros(N,1), zeros(N,1)],...
            [zeros(N,1),   cos(inc),  -sin(inc)],...
            [zeros(N,1),   sin(inc),   cos(inc)]);
QIC = permute(QIC, [3,2,1]);

end

