function status = compile_covertree_code()
%compile_covertree_code
%
%    This compiles the mex file for cover tree.  This is used for
%    cover tree approximate nearest neighbor search
%
% Copyright (c) 04-05-2017,  Zhouye Chen

Static_Dir = mfilename('fullpath');
Static_Dir = fileparts(Static_Dir);

Mesh_CPP    = fullfile(Static_Dir, 'src_code_covertree', 'mexCovertree.cpp');
Mesh_MEXDir = Static_Dir;
Mesh_MEX    = 'mexCovertree_CPP';
FLAG = 'CFLAGS=''\$CFLAGS -fopenmp'' LDFLAGS=''\$LDFLAGS -fopenmp''';% COPTIMFLAGS=''\$COPTIMFLAGS -fopenMP -02'' LDOPTIMFLAGS=''\$LDOPTIMFLAGS -fopenmp -02'' DEFINES=''\$DEFINES -fopenmp''';

disp('=======> Compile ''Cover tree''...');
mex -v CFLAGS='\$CFLAGS -fopenmp' LDFLAGS='\$LDFLAGS -fopenmp' ...
    COPTIMFLAGS='\$COPTIMFLAGS -fopenMP'...
    LDOPTIMFLAGS='\$LDOPTIMFLAGS -fopenmp'...
    DEFINES='\$DEFINES -fopenmp'...
    -largeArrayDims src_code_covertree/mexCovertree.cpp -output mexCovertree_CPP


% % % status = feval(@mex, '-v', '-largeArrayDims', Mesh_CPP, '-outdir', Mesh_MEXDir, '-output', Mesh_MEX);
% % % using omp
% % status = feval(@mex, '-v', FLAG, '-largeArrayDims', Mesh_CPP, '-outdir', Mesh_MEXDir, '-output', Mesh_MEX);
% % % use this for debug
% % % status = feval(@mex, '-g', '-v', '-largeArrayDims', Mesh_CPP,FLAG, '-outdir', Mesh_MEXDir, '-output', Mesh_MEX);
end