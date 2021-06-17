%% AD 02-02-2021
% Script to call AR_mypolyhedral.m function
% Requires use of AR's matlab scripts for functionality

%% Read in the .obj file

target='162p'

% root_dir: folder containing current work & pictures folder
root_dir='~/Documents/Year1/shape_modelling/162p/'
directory=[root_dir 'dist_cal_pole_scan/run_0_-90_2021_05_04/' ]

% Change to select which model is to be used - this part will be appended
% to file names!

% lam_beta:
modstr='300_40'
% object file:
objfilename = [directory target '_' modstr '.obj']
parfilename = [directory 'output_convex_pars_' modstr ]

%% Read-in the convex shape model

% Read vertices and facets:
[V,F]=read_vertices_and_faces_from_obj_file(objfilename);

% Calculate facet normals and facet areas
[FN,FNA]=AR_calcFN(V,F);

%% Set parameters rho and reff from literature:

rho = 4000 
reff=7 %km

%% Call AR_mypolyhedral function:

[results] = AR_mypolyhedral(V,F,rho,reff);

% print axial ratios a/b, b/c:

axial_ab = (results.a)/(results.b)
axial_bc = results.b/results.c
axial_ac = results.a/results.c
flattening_cb = results.c/results.b