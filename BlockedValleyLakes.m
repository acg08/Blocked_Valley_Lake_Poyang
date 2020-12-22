clear;
clc;
tstart=tic;

%%%%%%选择参数%%%%%%
%ep: small number;
%timeyear: coefficient convert year to second;
%g: gravitaional coefficient;
%rho: water density;
%R: submerged specific gravity;
%L: channel length;
%inis: initial slope;
%Duration: duration of calculation;
%If: flood intermittency;
%Q: flow discharge (m3/s);
%B: channel width;
%qw: flow discharge per unit width (m2/s);
%psi: ratio of tributrary and main stem bankfull discharge;
%Cz: dimensionless Chezy coefficient;
%Cf:dimensionless resistance coefficient;
%D: grain size of every group;
%p: bed porosity;
%au: coefficient related to Exner discretization;
%b: rate of base level rise (mm/year);
%afEH:adjustment factor of EH relation;
ep=1e-6;
If=0.73;
timeyear=365.25*24*3600;
timeyear=timeyear*If;
g=9.81;
rho=1000;
R=1.65;
L=1000000*4;
inis=0.00002115;
Duration=9000;
Nt=6;
%Main stem bankfull character
Q=38500;
B=1250;
qw=Q/B;
psi=0.109;
%Tributary bankfull character
Q=38500*psi;
B=B*psi^0.669;
qw=qw*0.48;
Cz=22.22;
Cf=Cz^(-2);
D=0.17/1000;
p=0.4;
au=1;
b=6;
afEH=1.7;

%%%%%%时空步长%%%%%%
%dt: time step (year);
%dx: cell size;
%n: node number;
dt=1;
dx=20000;
n=L/dx+1;

%%%%%%初始条件%%%%%%
%zb: bed elevation;
%s: bed slope;
x=linspace(0,L,n)';
zb=inis*(L-x);
zbini=zb;
s=[(zb(1)-zb(2))/dx;(zb(1:n-2)-zb(3:n))/2/dx;(zb(n-1)-zb(n))/dx];
%h: water depth;
%hini: initial water depth;
%u: flow velocity;
%taub: bed shear stress;
%shi: shields number;
%qstar: dimensionless bedload transport rate;
%qs: sediment transport rate per unit width;
%qsf: volume sediment supply per unit width (m2/s);
%Bl: incipient blocking parameter;
h=(Cf*qw.^2./g./s).^(1/3);
hini=h(n);
u=qw./h;
taub=rho.*Cf.*u.^2;
shi=taub./rho./R./g./D;
qstar=0.05./Cf.*shi.^2.5;
qstar=afEH*qstar;
qs=qstar.*(R*g*D)^0.5*D;
qsf=qs(1);
Gtf=qsf*B*3600*24*365.25*2650/1e9*If;
Bl=(1-p)*b/1000/365.25/24/3600*hini/qsf/inis/If;
t=0;
i=0;

%%%%%%定义变量%%%%%%
%Fr: Froude number;
%z: water surface;
Fr=zeros(n,1);
z=zeros(n,1);
%qsback: qs at the upstream node;
%qsit: qs at the node;
%qsfrnt: qs at the downstream node;
%qsdif: spatial difference of qs;
%dzb: change of bed evolution;
qsback=zeros(n,1);
qsit=zeros(n,1);
qsfrnt=zeros(n,1);
qsdif=zeros(n,1);
dzb=zeros(n,1);
%zbt: bed elevation at different time;
%st: bed slope at different time;
%qst: sediment transport rate per unit width at different time;
%ht: water depth at different time;
zbt=zeros(n,Nt);
st=zeros(n,Nt);
qst=zeros(n, Nt);
ht=zeros(n,Nt);
it=1;

%%%%%%Time Processing%%%%%%

while t<Duration
    %%%%%%Flow Hydraulics：Backwater Equation%%%%%%
    z(n)=hini+b*t/1000;
    h(n)=z(n)-zb(n);
    for ih=1:n-1
        hp=h(n-ih+1)-((zb(n-ih)-zb(n-ih+1))/dx-Cf*qw^2/g/h(n-ih+1)^3)/(1-qw^2/g/h(n-ih+1)^3)*dx;
        h(n-ih)=h(n-ih+1)-0.5*(((zb(n-ih)-zb(n-ih+1))/dx-Cf*qw^2/g/h(n-ih+1)^3)/(1-qw^2/g/h(n-ih+1)^3)+((zb(n-ih)-zb(n-ih+1))/dx-Cf*qw^2/g/hp^3)/(1-qw^2/g/hp^3))*dx;
    end    
    u=qw./h;
    Fr=u./(g.*h).^0.5;
    
    %%%%%%计算床面变形%%%%%%
    %%%sediment transport: Engelund Hansen%%%
    taub=rho.*Cf.*u.^2;
    shi=taub./rho./R./g./D; 
    qstar=0.05./Cf.*shi.^2.5;
    qstar=afEH*qstar;
    qs=qstar.*(R*g*D)^0.5*D;
    
    %%%bed evolution
    qsback=[qsf;qs(1:n-1)];
    qsit=qs;
    qsfrnt=[qs(2:n);2*qs(n)-qs(n-1)];
    qsdif=au*(qsit-qsback)+(1-au)*(qsfrnt-qsit);
    dzb=-qsdif/dx*dt*timeyear/(1-p);
    zb=zb+dzb;
    
    %%%parameter update
    t=t+dt;
    i=i+1;
    s=[(zb(1)-zb(2))/dx;(zb(1:n-2)-zb(3:n))/2/dx;(zb(n-1)-zb(n))/dx];
    z=h+zb;
    %%%record information every several years
    if (mod(i,Duration/dt/Nt)==0)
        zbt(:,it)=zb;
        st(:,it)=s;
        qst(:,it)=qs;
        ht(:,it)=h;
        it=it+1;
    end
    
end

tend=toc(tstart);
tc=tend/60