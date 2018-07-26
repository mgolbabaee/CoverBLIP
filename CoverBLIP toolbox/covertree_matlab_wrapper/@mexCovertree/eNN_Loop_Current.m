function [Label, D_atom, dist, counter] = eNN_Loop_Current(obj,pts,epsilon,upperbound)
%enn_loop_current: epsilon nearest neighbors search with a loop
% input:
%                      pts: enquery points, dxN, with d points of dimension N
%                epsilon: error of eNN search
%        upperbound: current estimation
% output:
%                   label: label of the enn search results in the dictionary
%               D_atom:  the atom of enn search results (dxN)
%                     dist:  the distance between the enqury points and the dictionary atom
%               counter:  the counted number for the distance calculations

if (nargin < 4)
    upperbound = realmax('single');
end

[Label, D_atom, dist, counter] = mexCovertree_CPP('enn_loop_current', obj.cppHandle, pts, epsilon, upperbound);

end