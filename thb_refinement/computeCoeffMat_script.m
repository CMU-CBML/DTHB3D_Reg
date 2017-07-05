% COMPUTECOEFFMAT_SCRIPT   Generate MEX-function computeCoeffMat_mex from
%  computeCoeffMat.
% 
% Script generated from project 'computeCoeffMat.prj' on 05-Jul-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'computeCoeffMat'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[27  1]);
ARGS{1}{5} = coder.typeof(0,[27  1]);
ARGS{1}{6} = coder.typeof(0,[27 64  2]);
ARGS{1}{7} = coder.typeof(0,[1 27]);

%% Invoke MATLAB Coder.
codegen -config cfg computeCoeffMat -args ARGS{1}

