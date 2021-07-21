%% Script to read in shape models (.obj files)
% created by AD 07-07-2021
close all

% root directory
root_dir = '/Users/s1523386/Documents/year1/shape_modelling/67p/';
date='20210720'

% Object name params
target = '67P';
lam = 79;
bet = 42;
namecore = '_50kfacets';
quick_plot = 1; % set to 1 if you want script to plot the model

%obj_filename = strcat(target,'_',int2str(lam),'_', int2str(bet),namecore,'.obj');

obj_filename = 'cg_spc_shap8_006k_cart.wrl'
obj_file = [ root_dir, obj_filename ];

% Read vertices and facets:
[V,F]=read_vertices_and_faces_from_wrl_file(obj_file,4);

% Change this if you're using wrl (starts at 0) vs obj (starts at 1)

% Calculate facet normals and facet areas
[FN,FNA]=AR_calcFN_wrl(V,F);
F = F+1;
%% Plot the model (if req at top of script i.e. quick_plot=1)
if  quick_plot==1
    % trisurf(T,x,y,z) plots the 3-D triangular surface defined by the points in vectors x, y, and z, and a triangle connectivity matrix T.
    figure
    trisurf(F,V(:,1),V(:,2),V(:,3))
    
    %colormap(pault_colormap('rb'))   % Use a custom colormap
    colormap gray
   % Ensure correct proportions:
     axis equal
end

%% Artificial LCs for non-convex shape model

objfilename=obj_filename;

% Calculate shadowing geometry matrices  
shadows=ShadowingGeometry(V,F,FN);

% Set up the spin-state completely manually

%model.t0 = 2452709 % From Rozek PhD thesis
model.t0 = 2456849.89034
model.lam = 79;
model.bet = 42;
model.p= 12.4041 % SHAP8 val
model.nu= 0;
model.omega = 0.045;
model.gF = -0.41;
model.B0 = 1.97;
model.hwidth  = 0.026;
model.rough = 15

%% Read in lightcurve

directory=root_dir;
lcfile='LC_LSST_3AU_NU';
lcfilename= strcat(directory, target,'_', lcfile, '.dat' );

% For the Mikko-format lightcurve (no uncertainties) use NUflag=1,
% magflag=0 (like data downloaded from DAMIT)

% For a file with uncertainties use NUflag=0,
% magflag=0 (our default when dealing with spin-state)

% For a file without uncertainties and with magnitudes use NUflag=1,
% magflag=1 (for example PH5 data for PH700)

NUflag=1;
magflag=0;

LCS_all = MikkoRead_v2(lcfilename,NUflag,magflag);

% How many lightcurves are read in?
lcmax=max(LCS_all(:,9))

% What is the peak-to-peak amplitude of each observed lightcurve?

 for lcno=1:lcmax
     lcrows=[];
     lcrows=find(LCS_all(:,9)==lcno);
     maxmag(lcno)=max(LCS_all(lcrows,10));
     minmag(lcno)=min(LCS_all(lcrows,10));
 end
 
%

amplmag=(maxmag-minmag)

%% Generate an example lightcurve 
%
lcid=6;
x=5; % Hapke scattering
ylims=1;

modstr = lcfile;

% With a non-convex shape
namecore=[target '_' lcfile '_' date '_Hapke']
lcres_one=Lightcurve_Generator_NonConv_v3(FN,FNA,shadows,LCS_all,model.t0,model.lam,model.bet,model.p,model.omega,model.B0,model.gF,model.hwidth,model.rough,x,model.nu,namecore,lcid,0,0.1,ylims, root_dir);