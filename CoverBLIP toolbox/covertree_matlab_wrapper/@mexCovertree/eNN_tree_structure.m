function tree_structure = eNN_tree_structure(obj)
% eNN_tree_structure: output tree structure in a matrix
% input:
%                      pts: enquery points, dxN, with d points of dimension N
% output:
%     tree_structure:  maximum level x number of atoms

tree_structure = mexCovertree_CPP('print_tree_structure', obj.cppHandle);

end