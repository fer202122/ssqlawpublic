%% norm_vecs: function description
function [norms] = norm_vecs(vectors, dim)
% vectors: N_vec * N_elem
norms = sqrt(sum(vectors.^2,dim));
end
