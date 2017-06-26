% STOREPIXELPHI_SCRIPT   Generate MEX-function storePixelPhi_mex from
%  storePixelPhi.
% 
% Script generated from project 'storePixelPhi.prj' on 24-Jun-2017.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'storePixelPhi'.
ARGS = cell(1,1);
ARGS{1} = cell(9,1);
ARGS{1}{1} = coder.typeof(int64(0));
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0,[Inf  3],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{6} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARGS{1}{7} = struct;
ARGS{1}{7}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{7}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{7}.chdElem = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{7}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{7}.node = coder.typeof(0);
ARGS{1}{7}.parElem = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7} = coder.typeof(ARGS{1}{7},[Inf  1],[1 0]);
ARGS{1}{8} = struct;
ARGS{1}{8}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{8} = coder.typeof(ARGS{1}{8},[Inf  1],[1 0]);
ARGS{1}{9} = struct;
ARGS{1}{9}.pU = coder.typeof(0);
ARGS{1}{9}.pV = coder.typeof(0);
ARGS{1}{9}.pW = coder.typeof(0);
ARGS{1}{9}.maxlevel = coder.typeof(0);
ARGS{1}{9}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.orderGauss = coder.typeof(0);
ARGS{1}{9}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.rho = coder.typeof(0,[3 1]);
ARGS{1}{9}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{9}.smallNumber = coder.typeof(0);
ARGS{1}{9}.lambda_1 = coder.typeof(0);
ARGS{1}{9}.lambda_2 = coder.typeof(0);
ARGS{1}{9} = coder.typeof(ARGS{1}{9});

%% Invoke MATLAB Coder.
codegen -config cfg storePixelPhi -args ARGS{1}

