% a function for 'readdata.m'
function[]=writefile(filename,row,Mr,Lr,a,e,T_ratio,gamma,F,precess,ps,U)
    [T_sat,narea,rarea]=calculate_Tsur_Ta(Mr,Lr,a,e,T_ratio,gamma,F,precess,ps,U);
    temp=[Mr,Lr,a,e,T_ratio,gamma,F,precess,ps,U,T_sat,narea,rarea];
    stemp=num2str(row);
    writematrix(temp,filename,'Range',strcat('A',stemp,':M',stemp));
end