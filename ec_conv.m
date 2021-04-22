%% Converting vectors from ecliptic to equatorial

function [ alph, dec ] = ec_conv(lam,bet);

e = 23.4393; % Obliquity of Earth's orbit in degrees

rot_mat = [ 1 0 0, 0 cosd(e) sind(e), 0 -sind(e) cosd(e) ];


alph = atand((sind(lam)*cosd(e)-tand(bet)*sind(e))/cosd(lam));
if (alph<0)
    alph = alph+360;
end
dec = asind((sind(bet)*cosd(e)+cosd(bet)*sind(e)*sind(lam)));


end
