%% AD 02-02-2021
% Script to call AR_mypolyhedral.m function
% Requires use of AR's matlab scripts for functionality
% 07-07-2021: Combined ^ and contents of basics tutorial so this script
% now:
% - Reads in .obj file
% - Calculates the axial ratios of the convex model (1)
% - Produces the publication-ready plot (2)
% - Generates the artificial lightcurves (3)

clear all;
clc;

% Run the first one, two or all three of these tasks by changing the ntask
% variable
ntask = 2;

target='162P'
date='20220428'

% root_dir: folder containing current work & pictures folder
% lc_generator saves images into root_dir/pictures/artifical_lcs
% Change to select which model is to be used - this part will be appended
% to file names!

root_dir='~/Documents/year1/shape_modelling/162p/'
directory = [ root_dir 'pole_scan/']%'pre_2022/dist_cal_pole_scan/run_0_-90_2021_06_01/' ] 

% %% Read in the lightcurve
%  
% lcdir = '~/Documents/year1/shape_modelling/162p/'
% lctarget = '162P'
% lcfile='LC_dist_cal_NU'
% lcfilename= strcat(lcdir, lctarget,'_', lcfile, '.dat' );
% 
% % Directory in which to save artificial lcs (2xN_LCs generated if lc=0)
% % alc_pics = [ root_dir 'pictures/' date '/artificial_lcs/' ]
% % 
% % if ~exist(alc_pics, 'dir')
% %     mkdir(alc_pics)
% % end

%% Read in the .obj file

% Change to select which model is to be used - this part will be appended
% to file names!

% lam_beta:1.5611
modstr='_100_-50';
% object file:
objfilename = [directory target modstr '.obj'];
%parfilename = [directory 'output_convex_pars_' modstr '_BS' ];

% If plotting nononvex wrl file: include n, the number of digits in the
% number of facets
% e.g. 3k facets: =4
n=4;

%% Read-in the convex shape model

% Read vertices and facets:

if (objfilename(length(objfilename)-2:length(objfilename))=='obj')
    [V,F]=read_vertices_and_faces_from_obj_file(objfilename)
else
    [V,F]=read_vertices_and_faces_from_wrl_file(objfilename,n);
    F = F+1;
end

% Calculate facet normals and facet areas
[FN,FNA]=AR_calcFN(V,F);


%% Set parameters rho and reff from literature:

rho = 4000;
reff=7; %km

%% Call AR_mypolyhedral function:
if ntask == 1 || ntask == 2 || ntask ==3
    [results] = AR_mypolyhedral(V,F,rho,reff);

    % print axial ratios a/b, b/c:

    axial_ab = (results.a)/(results.b)
    axial_bc = results.b/results.c
    axial_ac = results.a/results.c
    flattening_cb = results.c/results.b
end

%% AD: Quick tricolour plot of the convex surface


% % line 104: Change lcdir for Mikkoread
% % line 148 ish : Namecore for artificial lc changed to reflect nuflag=0 =>
% % take away _u if using lc 162P_LC_NU.dat
% 
% 
% % Plot the model, to see how it looks like. Use trisurf, e.g.: 
% figure
% trisurf(F,V(:,1),V(:,2),V(:,3))
% 
% % Use a custom colormap
% colormap(pault_colormap('rb'))
% 
% % Ensure correct proportions:
%  axis equal 
%  

%% Plot the model in the publication-ready format

 if ntask == 2 || ntask == 3
     % Set up some parameters
     diffstr=0.7; specularstr=0.1;
     alignx=1;


     %[ filename, plotfigure ] = plot_shape_publication(V, F, root_dir, target, modstr, date, diffstr, specularstr, alignx);
     plot_shape_publication(V, F, root_dir, target, modstr, date, diffstr, specularstr, alignx);

     % To plot with black background, uncomment the rest of the commands
     % below:
     
%      figure(plotfigure.Number)
%      %set(0,'DefaultAxesColor','black')
%      set(gca,'color','black')
%      %set(gcf, 'color','black')
%      
%      filename=[ filename '_dark']
%      
%      %set(0,'DefaultAxesColor','white')
%      
%      print(plotfigure,'-dpng','-r600',filename);
     %print(plotfigure,'-dpdf','-r600','-fillpage',filename);
     
     
 
 end

% %% Plot the model in the publication-ready format 
% % (with principal axis of inertia and X-axis aligned)
% 
 %modstr2=[ modstr '_alignx' ];
 %alignx=1;
 %plot_shape_publication(V, F, root_dir, target, modstr2, date, diffstr, specularstr, alignx);

% %% Plot a fancy shape projection
% % The shape will be visible from viewec direction
% % Illumination will come from lightvec direction
% 
% viewvec=[3 3 3];
% lightvec=[-1,1,1];
% 
% % An arbitrary colormap can be used
% colormapsel=pault_colormap('rb'); % see 
% %colormapsel='gray' % for smooth gray surface
% unitstr=''
% 
% plot_shape_fancy(V, F, target, modstr, date, lightvec, viewvec, colormapsel, unitstr);


 
%% Generating artificial lightcurves 

if ntask == 3
    % %Read-in the model parameters from convexinv output
     model=read_model_from_par_file(parfilename);

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

    amplmag=(maxmag-minmag)

    % Generate artificial lightcurves for all available lightcurves

     close all
     lcid=0
     x=4 % The Hapke and Lommel-Seelinger scattering isn't that different, so use LS
     yorp=0 % assume yorp=0
    % 
     namecore=[target '_' modstr '_' lcfile '_' date ]
    % 
     ylims=1;
    % 
    % 
    % % For convex model
     lcres=Lightcurve_Generator_Conv_v3(FN,FNA,LCS_all,model.t0,model.lam,model.bet,model.p,0,0,0,0,0,x,yorp,namecore,lcid,0,0,ylims,root_dir,alc_pics);
end

%close all

% %% Plotting plane-of-sky projections with correct illumination
% 
% % If, in the last section, lcid=0 was used, then here lcid should be set to
% % the ID of the lightcurve that we want plotted (if we plotted all 80
% % lightcurves and only wanted a plane-of-sky plot for lightcurve 3 then lcid = 3 )
%  lcid=1
%   
%  x=4
%  
%   
% % An arbitrary colormap can be used to highlight the illuminated regions
%  colormapsel=pault_colormap('Iri'); % lighter colours - illuminated, purple - in shadow
%  colormapsel='gray' % for smooth gray surface 
% %  
% %  
%  for n = (1:9:length(find(lcres.datapoints(:,7)==lcid))) 
%      plot_shape_POS_light(V,F,target, modstr, date, lcres, lcid, n, x,  colormapsel,model, root_dir)
%  end

