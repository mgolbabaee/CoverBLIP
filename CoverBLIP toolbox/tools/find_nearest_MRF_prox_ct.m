function [coeff, coeff_ind, D_atom, comp] = find_nearest_MRF_prox_ct(CT, X, epsilon,D, coeff_ind_old)
% Accelerated matched-filtering using cover tree's (1+epsilon)-ANNS
% approximate search
% % Inputs: 
%                 CT: cover tree structure
%                 X: query batch / mrf image 
%                 epsilon: cover tree's search inaccuracy parameter
%                 D: dictionary
%                 coeff_ind_old: current search indices (used for initializing the search)

% Outputs: 
%                 coeff: distances of nearest search result to the query batch
%                 coeff-ind: nearest atoms of D (indices) to the query batch X
%                 D: nearest atoms of D
%                 comp: computation cost (number of tree nodes traversed to find search result)
%
% (c) Mohammad Golbabaee, 2017

tmp_n=sqrt(sum(abs(X).^2,1));
%normalize data
X = bsxfun(@times, X, 1./(tmp_n+1e-8));

if sum(coeff_ind_old(:))>0
dmin = sqrt(sum(abs( (D.atoms(coeff_ind_old))- X  ).^2,1));
else
    dmin=2;
end

[coeff_ind, D_atom, coeff,comp] = CT.eNN_Loop_Current(transpose(X),epsilon,dmin); %new ct

if nargin>5
    %D_atom = D.atoms(coeff_ind);
    indd=find(coeff> dmin);
    if numel(indd)>0
        coeff(indd)= dmin(indd);
        coeff_ind(indd)=coeff_ind_old(indd);
        D_atom(:,indd) = (D.atoms(coeff_ind(indd))).';
    end
end

%% %find PD
coeff=sum(X.*conj(D_atom),1);
%coeff = (2 - coeff.^2)/2;

coeff = max(0, coeff.*tmp_n);
comp = comp*size(X,1);

end

