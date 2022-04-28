% Script that calls AD_UNCORRECT to output differences between mag zero and
% point on rotational lightcurve.

root_dir = '/Users/s1523386/Documents/year1/shape_modelling/67p/';
% Name of directory (inside root_dir) to save pictures in:
pics_file = 'pics_changeparams/'
% Declare if you want to plot everything you asked it to...
plot_artificial_lcs = 0;

% Object name params
target = '67P';
date='20211206_11';
%obj_filename = strcat(target,'_',int2str(lam),'_', int2str(bet),namecore,'.obj');

obj_filename = 'cg_spc_shap8_003k_cart.wrl'
obj_file = [ root_dir, obj_filename ];

% LC file name params
directory=root_dir;
lcfile='LC_LSST_ALL_NU_dummy';
lcfilename= strcat(directory, target,'_', lcfile, '.dat' );
NUflag=1;     % If you're using Mikkoformat w/o units (if dummy =1)
magflag=0;    % 0 if using intensities in Mikkofile
lcid=0;       % 0 for all LCs in file or number (e.g. 238) for that LC

%%%%%-----WRL OR OBJ FILE?-----%%%%%
% Read vertices and facets:

if obj_file(length(obj_file)-2:length(obj_file))=='wrl'
    [V,F]=read_vertices_and_faces_from_wrl_file(obj_file,4);
    % Change this if you're using wrl (starts at 0) vs obj (starts at 1)
    % Calculate facet normals and facet areas
    [FN,FNA]=AR_calcFN_wrl(V,F);
    F = F+1;
end

%% Set up the spin-state completely manually
%model.t0 = 2452709 % From Rozek PhD thesis
model.t0 = 2459953
model.lam = 79;
model.bet = 42;
model.p= 11.1 % SHAP8 val: 12.4041h
% Hapke scattering factors: Fornasier 2015 649nm (SPG)
model.omega = 0.1;
model.gF = -0.41;
model.B0 = 1.97;
model.hwidth  = 0.026;
model.rough = 15

T0 = model.t0;
lambda = model.lam;
beta = model.bet;
P = model.p;
omega = model.omega;
B0 = model.B0;
gF = model.gF;
hwidth = model.hwidth;
rough = model.rough;

ylims=1;
x = 4; % (3,4,5 = Lam,LS,Hapke: LEAVE AS HAPKE IF PLOTTING, script plots all models)

%% %results = AD_UNCORRECT(root_dir,date,target,obj_filename,V,F,FN,FNA,model.t0,model.lam,model.bet,model.p,model.omega,model.gF,model.B0,model.hwidth,model.rough,lcfile,lcfilename,NUflag,magflag,lcid,ylims,x,plot_artificial_lcs,pics_file)

results = AD_UNCORRECT(root_dir,date,target,obj_filename,V,F,FN,FNA,T0,lambda,beta,P,omega,gF,B0,hwidth,rough,lcfile,lcfilename,NUflag,magflag,lcid,ylims,x,plot_artificial_lcs,pics_file)

%% Plotting plane-of-sky projections with correct illumination

% If, in the last section, lcid=0 was used, then here lcid should be set to
% the ID of the lightcurve that we want plotted (if we plotted all 80
% % lightcurves and only wanted a plane-of-sky plot for lightcurve 3 then lcid = 3 )


% An arbitrary colormap can be used to highlight the illuminated regions
 colormapsel=pault_colormap('Iri'); % lighter colours - illuminated, purple - in shadow
 %colormapsel='gray' % for smooth gray surface 

% modstr='testing_cam'
% proj_dir = [ root_dir, 'projections/' ]
% for n = (1:9:length(find(results.datapoints(:,7)==1))) 
%      %print n
%      plot_shape_POS_light(V,F,target, modstr, date, results, 1, n, x, colormapsel,model, root_dir, proj_dir)
% end
