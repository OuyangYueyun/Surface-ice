% calculate incident stellar flux of particular point
function[F]=caculate_flux(fi,slat,h,r,L)
    cos_z=cos(fi)*cos(slat)*cos(h)+sin(fi)*sin(slat);
    if cos_z<0
        F=0;
    else
        S=L/(4*pi*r^2);
        F=S*cos_z;
    end
end