function sub_earth_lat = AD_coords(lam_p, bet_p, asc_node, omega, inc, a_init, d_init)

%% This function calculates the sub-Earth latitude of a given object
% Inputs:
% lam_p, bet_p: ecliptic long. and lat of pole orientation
% asc_node, omega: orbital ascending node and argument of perihelion
% inc: orbital inclination wrt ecliptic
% a_init, d_init: RA & dec of object as seen from Earth (this is reversed
% in the function to get position of Earth as seen from object)

% lam_p = 352.9076;
% bet_p = 63.2821;
% %lam_o = 319.5581;
% %bet_o = 88.1503;
eps0 = 23.4393;
% 
% asc_node = 49.558;
% omega = 286.502;
% inc = 1.850;
%      
% %a = 317.68;   % Pole RA
% %d = 52.89     % Pole dec
% [a,d] = ec_conv(lam_p,bet_p);
% a_init = 8.335;  % Earth-based equ. coords of Mars
% d_init = 3.6;
%% 
a_q = 180+a_init;   % Coords of person on planet looking back at Earth
d_q = -d_init;

% Matrix for transforming Earth-based eq -> ecliptic
C_cq = [1 0 0; 0 cosd(eps0) sind(eps0); 0 -sind(eps0) cosd(eps0)];

% Defining params for x vector:
% Rotation pole
z_Qc = [ cosd(lam_p)*cosd(bet_p); sind(lam_p)*cosd(bet_p); sind(bet_p)];
% Orbit pole:
z_Cc = [sind(asc_node)*sind(inc); -cosd(asc_node)*sind(inc); cosd(inc)];
%z_Cc = [cosd(lam_o)*cosd(bet_o); sind(lam_o)*cosd(bet_o); sind(bet_o) ];

% Orbital obliquity:
eps = acosd(dot(z_Cc, z_Qc)/(norm(z_Cc)*norm(z_Qc)));

%x vector:
x_Cc = cross(z_Qc, z_Cc)/abs(sind(eps));
x_Qc = x_Cc;

%y vectors:
y_Cc = cross(z_Cc, x_Cc);
y_Qc = cross(z_Qc, x_Qc);

%planet based ecliptic coordinates -> Earth based ecliptic coordinates
C_Cc = [ x_Cc'; y_Cc'; z_Cc' ];
C_cC = [ x_Cc'; y_Cc'; z_Cc' ]';

% Perifocus of planet's orbit (in Earth-based ecliptic coordinates):
P_c = [ (cosd(asc_node)*cosd(omega))-(sind(asc_node)*sind(omega)*cosd(inc)); (cosd(asc_node)*sind(omega)*cosd(inc))+(sind(asc_node)*cosd(omega)); sind(omega)*sind(inc) ];

P_C = C_Cc*P_c;

% Planet-based ecliptic longitude of the perifocus:
weird_pi = atan2(P_C(2), P_C(1));
weird_pi = rad2deg(weird_pi);

C_Qc = [ x_Qc'; y_Qc'; z_Qc' ];
C_cQ = C_Qc';

%C_Qq: Transform Earth-based eq to planet-Y based eq

C_Qq = C_Qc*C_cq;


%% Rectangular coords of unit vector in direction a_q, d_q relative to Earth based equatorial system

r_q = [ cosd(a_q)*cosd(d_q); sind(a_q)*cosd(d_q); sind(d_q) ];

% Transform this to planet-based Eq system:
r_Q = C_Qq*r_q;

lat_Q = asind(r_Q(3));
sub_earth_lat = lat_Q

end