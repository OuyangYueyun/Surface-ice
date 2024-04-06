% function used in 'readdata.m'
function[]=writefile(filename,row,Mr,Lr,a,e,T_ratio,gamma,F,precess,nyear)
    [narea,rarea]=caculate_Tsur(Mr,Lr,a,e,T_ratio,gamma,F,precess,nyear);
    temp=[Mr,Lr,a,e,T_ratio,gamma,F,precess,narea,rarea];
    stemp=num2str(row);
    writematrix(temp,filename,'Range',strcat('A',stemp,':J',stemp));
end