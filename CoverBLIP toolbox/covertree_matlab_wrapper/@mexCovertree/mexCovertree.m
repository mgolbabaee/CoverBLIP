%   MATLAB class wrapper for a Cover Tree C++ class
%
%   obj = mexCovertree(dict);
%
%   dict = MxN matrix representing the dictionary with M atoms

classdef mexCovertree < handle % this must be a sub-class of handle
    properties (SetAccess = private, Hidden = true)
        cppHandle         % Handle to the C++ class instance
    end
    properties (SetAccess = private)
        dict
    end
    methods
        % Constructor
        function obj = mexCovertree(varargin)
            
            if or((nargin < 1),(nargin > 1))
                disp('Requires 1 argument!');
                disp('The M-atom dictionary with N dimensions (MxN matrix).');
                error('Check the argument!');
            end
            
            obj.dict = varargin{1};
            % now create an instance of the C++ class
            % note: this will build the tree from the given dictionary
            obj.cppHandle = mexCovertree_CPP('new',obj.dict);
        end
        % Destructor - destroy the C++ class instance
        function delete(obj)
            mexCovertree_CPP('delete', obj.cppHandle);
        end
    end
end

% END %