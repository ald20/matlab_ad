function [ I, phi ] = AD_orbital_els(e,omega,asc_node,inc,q,lam,bet);
% Script to determine orbital obliquity (I) and argument of the sub-solar meridian at perihelion (phi)

%% Ecliptic coords orbit pole:
lam_n = 270+asc_node;
bet_n = 90-inc;

%% Ecliptic coords perihelion
tan_lam_p_long_asc = cosd(inc)*tand(omega);
lam_p = atand(tan_lam_p_long_asc)+asc_node;
if lam_p < 0
    lam_p = lam_p+360;
end

tan_bet_p = sind(lam_p-asc_node)*tand(inc);
bet_p = atand(tan_bet_p);

%% Ecliptic coords nu=90 degrees
tan_lam_v_long_asc = cosd(inc)*tand(omega+90);
lam_v = atand(tan_lam_v_long_asc)+asc_node;
if lam_v < 0
    lam_v = lam_v+360;
end

tan_bet_v = sind(lam_v-asc_node)*tand(inc);
bet_v = atand(tan_bet_v);

%% Angle of obliquity:

% Rotation axis omega_hat (om)
om = [ cosd(bet)*cosd(lam), cosd(bet)*sind(lam), sind(bet) ]';
om_norm = norm(om);

% Orbital momentum vector n_hat (n):
n = [ cosd(bet_n)*cosd(lam_n), cosd(bet_n)*sind(lam_n), sind(bet_n) ]';
n_norm = norm(n);

% Obliquity of small-body orbit (epsilon) (eps)
cos_eps = dot(n, om)/(n_norm*om_norm);
obliquity = acosd(cos_eps);
I = obliquity

%% Calculating argument of subsolar meridian at perihelion:

 cl = cosd(lam);     sl = sind(lam);
 cb = cosd(bet);     sb = sind(bet);
clp = cosd(lam_p);  slp = sind(lam_p);
cbp = cosd(bet_p);  sbp = sind(bet_p);
clv = cosd(lam_v);  slv = sind(lam_v);
cbv = cosd(bet_v); sbv = sind(bet_v);
cln = cosd(lam_n);  sln = sind(lam_n);
cbn = cosd(bet_n);  sbn = sind(bet_n);

w_vec = [ cb*cl, cb*sl, sb ]';
p_vec = [ cbp*clp, cbp*slp, sbp ]';
q_vec = [ cbv*clv, cbv*slv, sbv ]';
r_vec = [ cbn*cln, cbn*sln, sbn ]';

a = w_vec(1); b = w_vec(2); c = w_vec(3);
d = p_vec(1); e = p_vec(2); f = p_vec(3);
g = q_vec(1); h = q_vec(2); i = q_vec(3);
j = r_vec(1); k = r_vec(2); l = r_vec(3);

LHS = (a*h)-(b*g)-(h*j*cosd(I))+(k*g*cosd(I));
RHS = (e*g*sind(I))-(d*h*sind(I));

sin_phi = (LHS/RHS);

syms x;
phi_sols = solve(sind(x) == sin_phi, x);

sympref('FloatingPointOutput',true);
phi_sols;

phi = phi_sols


