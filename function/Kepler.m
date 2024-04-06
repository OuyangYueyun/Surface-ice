% calculate the velocity of revolution
function[E,r]=Kepler(J,a,e,k)
    r=caculate_r(a,e,k);
    E=J/(r*r);
end