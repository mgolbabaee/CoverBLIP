function [Label, D_atom, dist, counter] = eNN_Stop_Level(obj,pts,epsilon,stop_level)
%eNN_Stop_Level: epsilon nearest neighbors search and stop at a certain level
% input:
%                      pts: enquery points, dxN, with d points of dimension N
%                epsilon: error of eNN search
%            stop_level: stop searching level
% output:
%                   label: label of the enn search results in the dictionary
%               D_atom:  the atom of enn search results (dxN)
%                     dist:  the distance between the enqury points and the dictionary atom
%               counter:  the counted number for the distance calculations

if (nargin < 4)
    upperbound = realmax('single');
end

[Label, D_atom, dist, counter] = mexCovertree_CPP('enn_stop_level', obj.cppHandle, pts, epsilon, stop_level);

end