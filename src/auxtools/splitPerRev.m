function splitVar = splitPerRev(traj, fieldName)
%SPLITPERREV Splits fieldName to individual orbital revolutions
revs = traj.revs0;
wholeVar = traj.(fieldName);
splitVar = cell(revs,1);
idxStart = 1;
for ii=1:revs
    % Find end index of the current rev
    [~,idxEnd]=min(abs(traj.timerevs-ii));
    if isempty(idxEnd)
        error('Revolution not detected');
    elseif idxStart>idxEnd
        error('Start index larger than end index');
    end
    % Split
    splitVar{ii} = wholeVar(idxStart:idxEnd,:);
    % Update end index
    idxStart = idxEnd;
    
end

end