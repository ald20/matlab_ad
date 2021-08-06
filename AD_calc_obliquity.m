% Quick script to calculate obliquity and argument of sub solar meridian at perihelion of 67P
% Based on calculations from Rozek 2017 PhD thesis appendix B.2 & Sekanina
% 1981

%% Defining variables
% OSIRIS team rotation pole measurements:
% equatorial
%alpha = 69.3;
%dec = 64.1;

% ecliptic:
%lam = 78.1;
lam = 300;
bet = 38;
%bet = 41.5;

%% Orbital params
%(epoch 2015)
% e = 0.5404;
% omega = 308.58;
% asc_node = 214.53;
% inc = 6.3;
% q = 1.58; 


% % Orbital parameters (epoch 2021-03-29) (MPC)
%        e = 0.650017;  % eccentricity
%    omega = 22.096333; % Argument of periapsis (degrees)
% asc_node = 36.359973; % Longutude of ascending node (degrees)
%      inc = 3.871252;  % orbital inclination (degrees)
%        q = 1.210474;  % (AU)
       
%    % (epoch 2015-01-01) (Horizons)
%        e = 0.6409614;
%    omega = 12.78788;
% asc_node = 50.140206;
%      inc = 7.040258;
%        q = 1.243253;

% Orbital parameters 162P (Horizons) epoch 2010-Feb-27

e = 0.596083;
omega = 356.305763;
asc_node = 31.240098;
inc = 27.816583;
q = 1.233070;
   
[I, phi ] = AD_orbital_els(e, omega, asc_node, inc, q, lam, bet)