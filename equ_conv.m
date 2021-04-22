%% Converting vectors from equatorial to ecliptic

function [lam,bet] = equ_conv(alpha,dec);

e = 23.4; % Obliquity of Earth's orbit in degrees

rot_mat = [ 1 0 0, 0 cosd(e) sind(e), 0 -sind(e) cosd(e) ];


lam = atand((sind(alpha)*cosd(e)+tand(dec)*sind(e))/cosd(alpha));
bet = asind((sind(dec)*cosd(e)-cosd(dec)*sind(e)*sind(alpha)));

%coords_eq.lambda = lam;
%coords_eq.beta = bet;

end
