function [coeff, coeff_ind, D_atom, comp1] = find_nearest_MRF_mat(D, X, N, M)
% Exact matched-filtering with brute-force searches
% % Inputs: 
%                 X: query batch / mrf image 
%                 epsilon: cover tree's search inaccuracy parameter
%                 D: dictionary

% Outputs: 
%                 coeff: distances of nearest search result to the query batch
%                 coeff-ind: nearest atoms of D (indices) to the query batch X
%                 D: nearest atoms of D
%                 comp1: computation cost 
%
% (c) Mohammad Golbabaee, 2017
%%

coeff = zeros(1, N*M);

for ind_l = 1:N
    
    [coeff(1+(ind_l-1)*M:ind_l*M), coeff_ind(1+(ind_l-1)*M:ind_l*M)] = ...
        max(real(D.times(X(:,1+(ind_l-1)*M:ind_l*M))), [], 1);
    
end
D_atom = D.atoms(coeff_ind);
coeff = max(0, coeff);

natoms=size(D.times(ones(size(X,1),1)),1);
comp1 = natoms*numel(X);
end