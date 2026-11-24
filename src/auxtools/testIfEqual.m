function testIfEqual(vec,vecRef,tolOpt)
    % Check for optional input tolerance
    if nargin > 2
      tol = tolOpt;
    else
      tol = eps;
    end

    % Expand vecRef in matrix vs vector comparisons
    [rows,cols] = size(vec);
    if rows > 1 && (isrow(vecRef) || isscalar(vecRef))
        vecRef = vecRef.*ones(rows,cols);
    end

    % Assert if difference is below tolerance
    try
        diff = vec - vecRef;
        assert(all(abs(diff) < tol,'all'))
    catch ME
        switch ME.identifier
            case 'MATLAB:assertion:failed'
                [~,idx] = max(diff,[],"all","linear");
                [row,~] = ind2sub([rows,cols],idx);
                diffstr = num2str(diff);
                errID = 'compare:different';
                msg = ['Diff (max val): ', diffstr(row,:)];
                ME = addCause(ME,MException(errID,msg));
                throw(ME)
            otherwise
                rethrow(ME)
        end
    end

end

