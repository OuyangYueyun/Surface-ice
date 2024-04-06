function[tt,S_ice,f_str,ngrid]=calculate_iceage(T_sur,Ta,T_sat,ps,g,lat,lon,n1,n2,n3,nyear)
    % use the last cycle
    start=floor(n3*(nyear-4)/nyear);
    % params
    dens_ice=920;
    M_H2O=18e-3;
    d=1e3;
    R=8.31447;
    NA=6.02214e23;
    S_ice=0;
    S_esc=4*pi;
    Ta=mean(Ta);
    b_ia=1.9*10^21*(Ta/300)^0.75;
    t=T_sat-273.15;
    p1=exp(43.494-6545.8/(t+278))/(t+868)^2;
    f_str=p1/ps;
    Fi=f_str*b_ia*(28-2)*1e-3*g/(R*Ta*NA);
    
    ngrid=0;
    area=caculate_rarea(lat,lon);       % grid area
    for i=1:n1
        for j=1:n2
            flag=1;
            for t=start:n3
                if T_sur(i,j,t)>T_sat
                    flag=0;
                    break;
                end
            end
            if flag==1
                ngrid=ngrid+1;
                S_ice=S_ice+area(i,j);
            end
        end
    end
    n_ice=dens_ice*S_ice*d/M_H2O;
    tt=n_ice/(S_esc*Fi*365*1e9*23.97*3600);     % duration of ice sheet coverage
end