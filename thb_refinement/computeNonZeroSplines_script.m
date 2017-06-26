% COMPUTENONZEROSPLINES_SCRIPT   Generate MEX-function computeNonZeroSplines_mex
%  from computeNonZeroSplines.
% 
% Script generated from project 'computeNonZeroSplines.prj' on 21-Jun-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'computeNonZeroSplines'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{2} = struct;
ARGS{1}{2}.pU = coder.typeof(0);
ARGS{1}{2}.pV = coder.typeof(0);
ARGS{1}{2}.pW = coder.typeof(0);
ARGS{1}{2}.maxlevel = coder.typeof(0);
ARGS{1}{2}.nelemU = coder.typeof(0);
ARGS{1}{2}.nelemV = coder.typeof(0);
ARGS{1}{2}.nelemW = coder.typeof(0);
ARGS{1}{2}.orderGauss = coder.typeof(0);
ARGS{1}{2}.kU = coder.typeof(0);
ARGS{1}{2}.kV = coder.typeof(0);
ARGS{1}{2}.kW = coder.typeof(0);
ARGS{1}{2}.nobU = coder.typeof(0);
ARGS{1}{2}.nobV = coder.typeof(0);
ARGS{1}{2}.nobW = coder.typeof(0);
ARGS{1}{2}.rho = coder.typeof(0,[3 1]);
ARGS{1}{2}.timestep = coder.typeof(0,[3 1]);
ARGS{1}{2}.smallNumber = coder.typeof(0);
ARGS{1}{2}.lambda_1 = coder.typeof(0);
ARGS{1}{2}.lambda_2 = coder.typeof(0);
ARGS{1}{2} = coder.typeof(ARGS{1}{2});
ARGS{1}{3} = struct;
ARGS{1}{3}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{3}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{3}.chdElem = coder.typeof(0);
ARGS{1}{3}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{3}.node = coder.typeof(0);
ARGS{1}{3}.parElem = coder.typeof(0);
ARGS{1}{3}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{3});
ARGS{1}{4} = struct;
ARGS{1}{4}.basis_ind = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{4}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4}.chdBasis = coder.typeof(0);
ARGS{1}{4}.coeffBasis = coder.typeof(0);
ARGS{1}{4}.suppCell = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{4}.flag_trunc = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4}.flag_ref = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4}.flag_bdy = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4}.actB = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{4} = coder.typeof(ARGS{1}{4});

%% Invoke MATLAB Coder.
cd('C:\Users\cosmin\Dropbox\matlab\aishwarya\DTHB3D_ver4\THBrefinement');
codegen -config cfg computeNonZeroSplines -args ARGS{1}

