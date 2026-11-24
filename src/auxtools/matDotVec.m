function vecOut = matDotVec(Q,vecIn)
%MATDOTVEC Computes the matrix product between a time-series
% of rotation matrices and a vector
    N = size(vecIn,1);
    if isrow(vecIn)
        vecInPages = vecIn';
        vecOutPages = pagemtimes(Q,vecInPages);
        vecOut = vecOutPages';

    else
        vecInPages = reshape(vecIn', [3,1,N]);
        vecOutPages = pagemtimes(Q,vecInPages);
        vecOut = reshape(vecOutPages, [3,N])';

    end

end