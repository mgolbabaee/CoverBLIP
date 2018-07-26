function [label, D_atom, dist, counter, index] = eNN_BSearch(obj,pts,epsilon,batch_size)
% eNN_BSearch: epsilon nearest neighbors batch search
% input:
%                      pts: enquery points, dxN, with d points of dimension N
%                epsilon: error of eNN search
%           batch_size: batch search size, dividing d enquery points into batches
% output:
%                   label: label of the enn search results in the dictionary
%               D_atom:  the atom of enn search results (dxN)
%                     dist:  the distance between the enqury points and the dictionary atom
%               counter:  the counted number for the distance calculations
%                  index:   the original index of the input enquery points

[label, D_atom, dist, counter, index] = mexCovertree_CPP('enn_search_batch', obj.cppHandle, pts, epsilon, batch_size);

end