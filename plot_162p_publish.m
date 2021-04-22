%% Basic handling of models and lightcurves
% Requires AR matlab scripts
% Clean the workspace (use command+enter to evaluate current cell)

close all
clear all

%% Read in the shape and calculate shadowing (if you're going to use a non-convex shape model)

target='162P'
date='20201211'

% CHANGEMEEEEEEEEE

directory='polescan_fixed_32_864052/run_150_-90_2020_12_09_19_16_42/' 
modstr='230_-70'
objfilename = [directory target '_' modstr '.obj']
% parfilename = [directory 'output_parameters_run_' modstr '_AR']

% Read vertices and facets:
[V,F]=read_vertices_and_faces_from_obj_file(objfilename);
% Calculate facet normals and facet areas
[FN,FNA]=AR_calcFN(V,F);

% Setup some parameters
diffstr=0.7;
specularstr=0.1;
alignx=0;
plot_shape_publication(V, F, target, modstr, date, diffstr, specularstr, alignx);

close all

% Get misalignment stuff

rho=370 % density 67P
Reff=2010 % eff. radius 67P

[results] = AR_mypolyhedral(V,F,rho,Reff);

results
