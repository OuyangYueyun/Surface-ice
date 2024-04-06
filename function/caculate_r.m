% calculate the current distance from planet to star
function[r]=caculate_r(a,e,k)
    r=a*(1-e^2)/(1+e*cos(k));
end