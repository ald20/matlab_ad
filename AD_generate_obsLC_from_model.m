% Script that calls AD_UNCORRECT to output differences between mag zero and
% point on rotational lightcurve.

root_dir = '/Users/s1523386/Documents/year1/shape_modelling/67p/';
date='20210924_3k';
% Object name params
target = '67P';

%obj_filename = strcat(target,'_',int2str(lam),'_', int2str(bet),namecore,'.obj');

obj_filename = 'cg_spc_shap8_003k_cart.wrl'
obj_file = [ root_dir, obj_filename ];

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
model.p= 12.4041 % SHAP8 val
% Hapke scattering factors: Fornasier 2015 649nm (SPG)
model.omega = 0.045;
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

lcfile='LC_LSST_ALL_NU_dummy';
directory=root_dir;
lcfilename= strcat(directory, target,'_', lcfile, '.dat' );
NUflag=1;
magflag=0;
lcid=238; % 0 for all poss or number for that LC
ylims=1;
x = 3; % Lommel-Seeliger scattering (3,4,5 = Lam,LS,Hapke)

% Declare if you want this to plot everything you asked it to...
plot_artificial_lcs = 1;

results = AD_UNCORRECT(root_dir,date,target,obj_filename,V,F,FN,FNA,model.t0,model.lam,model.bet,model.p,model.omega,model.gF,model.B0,model.hwidth,model.rough,lcfile,lcfilename,NUflag,magflag,lcid,ylims,x,plot_artificial_lcs)