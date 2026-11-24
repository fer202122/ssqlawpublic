function output = linearAngle(input)
%LINEARANGLE Linearly unwraps an angle array
input = wrapTo2Pi(input);
idx = find(abs(diff(input))>deg2rad(355)) + 1;
for ii=1:length(idx)
    input(idx(ii):end) = input(idx(ii):end) + 2*pi;
end
output = input;

end