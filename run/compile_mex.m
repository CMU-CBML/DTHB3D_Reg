%script from automatically generating compiled *.mex files from the matlab
%*.m files using the Matlab Coder toolbox (codegen command). 


%adding paths
addpaths();


%% Create configuration object of class 'coder.MexCodeConfig'.
cfg = coder.config('mex');
cfg.GenerateReport = true;
%cfg.ReportPotentialDifferences = false;
cfg.GlobalDataSyncMethod = 'NoSync';

cd ./iteration_loop_func
%% Define argument types for entry-point 'computenewPoints'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{3} = struct;
ARGS{1}{3}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{3},[Inf  1],[1 0]);
ARGS{1}{4} = struct;
ARGS{1}{4}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{4} = coder.typeof(ARGS{1}{4},[Inf  1],[1 0]);
ARGS{1}{5} = struct;
ARGS{1}{5}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{5} = coder.typeof(ARGS{1}{5},[Inf  1],[1 0]);
ARGS{1}{6} = struct;
ARGS{1}{6}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{6} = coder.typeof(ARGS{1}{6},[Inf  1],[1 0]);
ARGS{1}{7} = coder.typeof(0);

%% Invoke MATLAB Coder.
codegen -config cfg computenewPoints -args ARGS{1}


%% Define argument types for entry-point 'compute_Integ_Domain'.
ARGS = cell(1,1);
ARGS{1} = cell(24,1);
ARGS{1}{1} = struct;
ARGS{1}{1}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS{1}{1},[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{3} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4  4],[1 0 0]);
ARGS{1}{5} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{6} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{7} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{8} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{9} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{10} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{11} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{12} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{13} = coder.typeof(single(0),[Inf  4  4],[1 0 0]);
ARGS{1}{14} = coder.typeof(0,[Inf  4],[1 0]);
ARGS{1}{15} = struct;
ARGS{1}{15}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{15} = coder.typeof(ARGS{1}{15},[Inf  1],[1 0]);
ARGS{1}{16} = struct;
ARGS{1}{16}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{16} = coder.typeof(ARGS{1}{16},[Inf  1],[1 0]);
ARGS{1}{17} = struct;
ARGS{1}{17}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{17} = coder.typeof(ARGS{1}{17},[Inf  1],[1 0]);
ARGS{1}{18} = struct;
ARGS{1}{18}.mat = coder.typeof(single(0),[Inf  4  4  4],[1 0 0 0]);
ARGS{1}{18} = coder.typeof(ARGS{1}{18},[Inf  1],[1 0]);
ARGS{1}{19} = coder.typeof(0);
ARGS{1}{20} = coder.typeof(0);
ARGS{1}{21} = coder.typeof(0,[4 1]);
ARGS{1}{22} = coder.typeof(0,[4 1]);
ARGS{1}{23} = coder.typeof(0,[4 1]);
ARGS{1}{24} = coder.typeof(single(0),[Inf  3],[1 0]);

%% Invoke MATLAB Coder.
codegen -config cfg compute_Integ_Domain -args ARGS{1}

%% Define argument types for entry-point 'tripleIterLoop'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0,[1 3]);
ARGS{1}{2} = struct;
ARGS{1}{2}.active_cell = coder.typeof(0);
ARGS{1}{2}.phi = coder.typeof(single(0),[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARGS{1}{3} = struct;
ARGS{1}{3}.nzsplines = coder.typeof(int64(0),[Inf  1],[1 0]);
ARGS{1}{3} = coder.typeof(ARGS{1}{3},[Inf  1],[1 0]);
ARGS{1}{4} = coder.typeof(0,[Inf  4],[1 0]);

%% Invoke MATLAB Coder.
codegen -config cfg tripleIterLoop -args ARGS{1}

%% Define argument types for entry-point 'storePixelPhi'.
cd ../../store_phi_functions

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

%% Define argument types for entry-point 'GaussPhi'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0,[Inf  2],[1 0]);
ARGS{1}{2} = struct;
ARGS{1}{2}.knot_ind = coder.typeof(0,[Inf  3  2],[1 0 0]);
ARGS{1}{2}.flag_active = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.IEN = coder.typeof(0,[Inf  27],[1 0]);
ARGS{1}{2}.chdElem = coder.typeof(0,[Inf  8],[1 0]);
ARGS{1}{2}.cell_centre = coder.typeof(0,[Inf  3],[1 0]);
ARGS{1}{2}.node = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.parElem = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2}.actE = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{2} = coder.typeof(ARGS{1}{2},[Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{3} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARG = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof({ARG}, [Inf  1],[1 0]);
ARGS{1}{6} = struct;
ARGS{1}{6}.mat = coder.typeof(single(0),[Inf  27],[1 0]);
ARGS{1}{6} = coder.typeof(ARGS{1}{6},[Inf  1],[1 0]);
ARGS{1}{7} = struct;
ARGS{1}{7}.pU = coder.typeof(0);
ARGS{1}{7}.pV = coder.typeof(0);
ARGS{1}{7}.pW = coder.typeof(0);
ARGS{1}{7}.maxlevel = coder.typeof(0);
ARGS{1}{7}.nelemU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nelemV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nelemW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.orderGauss = coder.typeof(0);
ARGS{1}{7}.kU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.kV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.kW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobU = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobV = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.nobW = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.rho = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.timestep = coder.typeof(0,[Inf  1],[1 0]);
ARGS{1}{7}.smallNumber = coder.typeof(0);
ARGS{1}{7}.lambda_1 = coder.typeof(0);
ARGS{1}{7}.lambda_2 = coder.typeof(0);
ARGS{1}{7} = coder.typeof(ARGS{1}{7});

%% Invoke MATLAB Coder.
codegen -config cfg GaussPhi -args ARGS{1}

cd ../bspline-util

%% Define argument types for entry-point 'FindSpan'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
codegen -config cfg FindSpan -args ARGS{1}

cd ../thb_refinement
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



cd ../run
