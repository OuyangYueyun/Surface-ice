% This code calculate the surface temperature of Mercury
clear;
addpath(genpath('..\function'))
% calculate the orbit (use radians!!!)
% set parameters
T_orb=87.9691*86400;    % revolution period
T_rot=58.646*86400;     % rotation period
M_day=86400/(86400/T_rot-86400/T_orb);        % length of Mercury day
a=57909050e3;    % semi-major axis
e=0.205630;    % eccentricity
dt=T_orb/(200);   % time step
b=a*sqrt(1-e^2);    % semi-minor axis
J=2*pi*a*b/T_orb;   % angular momentum
w=2*pi/(58.646*86400);   % spin velocity
gama=deg2rad(0.034);     % obliquity
precess=0;   % precession angle

% calculate time sequence of seasonal angle and longitude of substellar point
kappa=0;
slon=pi;
r=caculate_r(a,e,kappa);
% Euler's forward 
[E,~]=Kepler(J,a,e,kappa);
kappa(2)=kappa(1)+dt*E;
slon(2)=slon(1)-dt*(w-E);
while slon(2)<0
    slon(2)=slon(2)+2*pi;
end
% Central Difference
nyear=6;
while kappa(end)<kappa(1)+2*pi*nyear
    n=size(kappa,2);
    [E,r(n)]=Kepler(J,a,e,kappa(n));
    kappa(n+1)=kappa(n-1)+2*dt*E;
    slon(n+1)=slon(n-1)-2*dt*(w-E);
    while slon(n+1)<0
        slon(n+1)=slon(n+1)+2*pi;
    end
end
r(n+1)=caculate_r(a,e,kappa(end));

% calculate incident stellar flux of each point
% set parameters
n1=11;     % latitude grid point
n2=21;     % longitude grid point
n3=n+1;     % number of time steps
L=3.828e26;     % Luminosity of Sun
slat=asin(cos(kappa-precess).*sin(gama));     % latitude of substellar plot

lat=linspace(-pi/2,pi/2,n1);
lon=linspace(0,2*pi,n2);
flux=zeros(n1,n2,n3);        % (lat, lon, time)
for i=1:n1
    for j=1:n2
        for k=1:n3
            flux(i,j,k)=caculate_flux(lat(i),slat(k),lon(j)-slon(k),r(k),L);
        end
    end
end
% clearvars -except flux n1 n2 n3 dt slon

% calculate temperature of each soil layer
% params
absorp=1-0.088;
Kt=0.01;
rcp=1e6;     % lunar regolith
F=0;        % Internal heat flux
nz=35;
z=zeros(nz+1,1);
z0=sqrt(Kt*dt/rcp);
for i=1:(nz+1)
    z(i)=z0*(exp((i-1)/5)-1)/(exp(3)-1);
end
dz=zeros(nz,1);
for i=1:nz
    if i<nz
        dz(i)=(z(i+2)-z(i))/2;
    else
        dz(nz)=z(i+1)-z(i);
    end
end
T=zeros(n1,n2,n3,nz+1);
T(:,:,1,:)=100;

% solution
m=rcp/dt;
% caculate coef matrix A
A=zeros(nz+1);  
c2=-Kt/(z(2)-z(1));
A(1,2)=c2;
for j=2:nz
    c1=Kt/(dz(j-1)*(z(j)-z(j-1)));
    c3=Kt/(dz(j-1)*(z(j+1)-z(j)));
    c2=-(m+c1+c3);
    A(j,j-1:j+1)=[c1,c2,c3];
end
c1=Kt/(dz(nz)*(z(nz+1)-z(nz)));
c2=-(m+c1);
A(nz+1,nz:nz+1)=[c1,c2];

for lat=1:n1
    for lon=1:n2
        for i=2:n3
            %time-related terms in A
            c1=OLR(T(lat,lon,i-1,1),1)+Kt/(z(2)-z(1));
            A(1,1)=c1;
            %caculate B
            rad_in=flux(lat,lon,i);
            b0=absorp*rad_in-OLR(T(lat,lon,i-1,1),0)+OLR(T(lat,lon,i-1,1),1)*T(lat,lon,i-1,1);
            temp=squeeze(T(lat,lon,i-1,2:nz));
            b=cat(1,b0,-m*temp,-m*squeeze(T(lat,lon,i-1,nz+1))-F/dz(nz));
            %solve AX=B
            T(lat,lon,i,:)= tridiagonal(A,b)';
        end
    end
end
T_sur=T(:,:,:,1);
% clearvars -except T_sur n1 n2 n3

nd=M_day/dt;
%% line plot: surface temperature of particular point
t=(1:n3);
plot(t/nd-1.7,squeeze(T_sur(6,6,t)),'LineWidth',2)
grid on
xlim([0 1])
xlabel('Mercury day')
ylabel('Surface temperature / K')
width=550;
height=300;
left=100;
bottom=100;
set(gcf,'position',[left,bottom,width,height])

%% GIF: evolution of surface temperature distribution
clearvars -except T_sur n1 n2 n3
lat=linspace(-90,90,n1);
lon=linspace(-180,180,n2);
step=2;
fig = figure;
writerObj=VideoWriter('gif_Mercury.avi');
open(writerObj);
m_proj('robinson','lat',[-90 90],'lon',[-180 180]);
for i = 20:step:n3
    temp=T_sur(:,:,i);
    colormap([m_colmap('diverging')]);  
    [~,h]=m_contourf(lon,lat,temp,7);
    set(h,'Color','none');
    clim([50 700]);
    colorbar;
    title('Simulation of Mercury''s surface temperature','FontSize',14)
    m_grid;
    xlabel('lat')
    ylabel('lon')
    pause(0.01);
    frame = getframe(fig);
    writeVideo(writerObj,frame);
end
close(writerObj);