% This code calculate surface ice coverage under different conditions
% according to 'params.xlsx'.

clear;
filename='params.xlsx';
data=readmatrix(filename,'Range','A2:I5');
% Star type
temp=find(isnan(data(:,1)));
ns=temp(1)-1;
M=data(1:ns,1);     % Stellar mass
L=data(1:ns,2);      % Stellar Luminosity
% semi-major axis
na=2;
a_earth=149598023e3;
a_mercury=57909050e3;
% eccentricity
temp=find(isnan(data(:,3)));
ne=temp(1)-1;
e=data(1:ne,3);      
% revolution period: rotation period
temp=find(isnan(data(:,4)));
nt=temp(1)-1;
T_ratio=data(1:nt,4);      
% obliquity
temp=find(isnan(data(:,5)));
ng=temp(1)-1;
gamma=deg2rad(data(1:ng,5));      
% internal heat flux
temp=find(isnan(data(:,6)));
nf=temp(1)-1;
F_internal=1e-3*data(1:nf,6);          
% surface pressure
temp=find(isnan(data(:,8)));
np=temp(1)-1;
ps=data(1:np,8);
% surface wind speed
temp=find(isnan(data(:,9)));
nu=temp(1)-1;
U=data(1:np,9);
N=[ns na nt ng nf np nu];

%%
% reference state
M0=M(1);
L0=L(1);
a0=a_earth*sqrt(L0);
e0=e(1);
Tr0=T_ratio(1);
gamma0=gamma(1);
F0=F_internal(1);
precess0=0;
ps0=ps(1);
U0=U(1);

filename='result.xls';
string={'M','L','a','e','T_ratio','gamma','F','precess',...
    'ps','U','T_sat','ngrid','area_ratio'};
writecell(string,filename,'Range','A1:M1','WriteMode','overwritesheet');
writefile(filename,2,M0,L0,a0,e0,Tr0,gamma0,F0,precess0,ps0,U0);
row=3;

% change star
for i=2:ns
    writefile(filename,row,M(i),L(i),a_earth*sqrt(L(i)), ...
        e0,Tr0,gamma0,F0,precess0,ps0,U0);
    row=row+1;
end
% change semi-major axis
writefile(filename,row,M0,L0,a_mercury*sqrt(L0),e0,Tr0,gamma0,F0,precess0,ps0,U0);
row=row+1;
% change revolution period: rotation period
for i=2:nt
    writefile(filename,row,M0,L0,a_earth*sqrt(L0),e(2),T_ratio(i),gamma0,F0,precess0,ps0,U0);
    row=row+1;
end
% change obliquity
for i=2:ng
    writefile(filename,row,M0,L0,a_earth*sqrt(L0),e0,Tr0,gamma(i),F0,precess0,ps0,U0);
    row=row+1;
end
% change internal heat flux
for i=2:nf
    writefile(filename,row,M0,L0,a_earth*sqrt(L0),e0,Tr0,gamma0,F_internal(i),precess0,ps0,U0);
    row=row+1;
end
% change surface pressure
for i=2:np
    writefile(filename,row,M0,L0,a_earth*sqrt(L0),e0,Tr0,gamma0,F0,precess0,ps(i),U0);
    row=row+1;
end
% change surface wind speed
for i=2:np
    writefile(filename,row,M0,L0,a_earth*sqrt(L0),e0,Tr0,gamma0,F0,precess0,ps0,U(i));
    row=row+1;
end
%% draw figure
filename='result.xls'; row=15; N=[3,2,3,3,2,3,3];
clearvars -except filename row N
data=readmatrix(filename,'Range',strcat('M2:M',num2str(row-1)));
nparam=size(N,2);
nd=max(N);
Y=NaN(nd,nparam);
X=1:nparam;
count=2;

figure
width=650;
height=400;
left=100;
bottom=100;

% scatter
for i=1:nparam
    Y(1,i)=data(1);
    Y(2:N(i),i)=data(count:count+N(i)-2);
    scatter(ones(size(Y(:,i)))*i,Y(:,i),100,'MarkerFaceColor','b',...
        'MarkerEdgeColor','k');
    hold on
    count=count+N(i)-1;
end
% dash line (ice coverage of reference state)
plot([0 10],[data(1),data(1)],'k--','LineWidth',2);
xlim([0.5 nparam+0.5]);
ylim([0 1])
box on
xticks(1:9);
xticklabels({'Star','F_{star}','e+T_{ratio}','\gamma','F_{internal}', ...
    'p_0','U'})
yticks(0:0.2:1);
ylabel('Surface ice coverage')
grid on
hold on
% deviation from reference state
ymax=max(Y);
yneg=range(Y,1);
ypos=zeros(size(ymax));
errorbar(X,ymax,yneg,ypos,'r','LineStyle','none','LineWidth',2,'CapSize',0);
set(gcf,'position',[left,bottom,width,height])
hold off
