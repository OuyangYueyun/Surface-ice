
function[h]=terminator(phi,delt)
    temp=-tan(phi)*tan(delt);
    if temp<-1
        h=pi;
    elseif temp>1
        h=0;
    else
        h=acos(temp);
    end
end