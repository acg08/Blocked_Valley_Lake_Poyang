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
%Phi: phi value of grain size;
%D: characteristic grain size of each size range;
%ng: number of sediment group;
%p: bed porosity;
%au: coefficient related to Exner discretization;
%blr: rate of base level rise (mm/year);
%a,b,c,d: parameters in SEH;
%Ls: thickness of storage layer;
%altr: coefficient related to interfacial exchange;
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
Q=38500; %main stem
B=1250;  %main stem
qw=Q/B;
qw=qw*0.48;
Cz=22.22;
Cf=Cz^(-2);
Phin=[-4,-3,-2,-1,0,1];
Phi=(Phin(1:5)+Phin(2:6))/2;
D=2.^Phi/1000;
ng=size(D,2);
p=0.4;
au=1;
blr=6;
a=0.037;
b=0.445;
c=0.876;
d=-0.348;
Ls=0.5;
altr=0.5;
afSEH=2.3;

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
%Fi: sediment fraction of bed surface;
%Dsg: geometric mean grain size of bed surface;
x=linspace(0,L,n)';
zb=inis*(L-x);
zbini=zb;
s=[(zb(1)-zb(2))/dx;(zb(1:n-2)-zb(3:n))/2/dx;(zb(n-1)-zb(n))/dx];
Fi=ones(n,1)*[13.8,76.1,7.7,1.2,1.2]/100;
Dsg=2.^(Fi*Phi')/1000;
%h: water depth;
%hini: initial water depth;
%u: flow velocity;
h=(Cf*qw.^2./g./s).^(1/3);
hini=h(n);
u=qw./h;
%taub: bed shear stress;
%ustar: shear velocity;
%shi: shields number for sediment of each size range;
%NA: coefficient in Naito's relation;
%NB: exponent in Naito's relation;
%Nstar: dimensionless sediment transport rate;
%qse: equilibrium sediment transport rate per unit width for each size range;
%qs: sediment transport rate per unit width;
%qsT: total volume sediment transport rate per unit width;
%Gta: ambient annual sediment load;
%qsf: sediment supply rate per unit width;
%qsfT: total sediment supply rate per unit width;
%Gtf: annual sediment supply;
%Bl: incipient blocking parameter;
taub=rho.*Cf.*u.^2;
ustar=u.*(Cf)^0.5;
shi=taub*ones(1,ng)./rho./R./g./(ones(n,1)*D);
NA=a.*((ones(n,1)*D)./(Dsg*ones(1,ng))).^b;
NB=c.*((ones(n,1)*D)./(Dsg*ones(1,ng))).^d;
Nstar=NA.*(shi).^NB;
qs=Nstar.*Fi.*((ustar.^3)*ones(1,ng))./R./g./Cf;
qs=qs*afSEH;
qsT=sum(qs,2);
Pstr=qs./(qsT*ones(1,ng));
Dg_load=2.^(Pstr*Phi')/1000;
Gta=qsT/1000*timeyear*rho*(R+1)*B/1e6;
qsf=qs(1,:);
qsfT=qsT(1);
Gtf=qsfT*B*3600*24*365.25*2650/1e9*If;
Bl=(1-p)*blr/1000/365.25/24/3600*hini/qsfT/inis/If;
%La: thickness of active layer;
%Psub: GSD of substrate sediment;
%Store: information of substrate stratigraphy;
%indup: number of sublayer at every node;
%Pup:proportion of sediment on uppermost sublayer;
%Lup: thickness of uppermost sublayer
La=ones(n,1)*0.2*hini;
Psub=(Pstr(1,:)+Fi(1,:))/2;
Store=cell(n,1);
indup=ceil((zb-La+10*Ls)./Ls);
for is=1:n
    Store{is,1}=zeros(300,ng);
    Store{is,1}(1:indup(is),:)=ones(indup(is),1)*Psub;
end
Pup=ones(n,1)*Psub;
Lup=zb+10*Ls-La-Ls*(indup-1);
Store_ini=Store;
t=0;
i=0;

%%%%%%定义变量%%%%%%
%Fr: Froude number;
%z: water surface;
Fr=zeros(n,1);
z=zeros(n,1);
%qsTback: qsT at the upstream node;
%qsTit: qsT at the node;
%qsTfrnt: qsT at the downstream node;
%qsTdif: spatial difference of qsT;
%dzb: change of bed elevation;
qsTback=zeros(n,1);
qsTit=zeros(n,1);
qsTfrnt=zeros(n,1);
qsTdif=zeros(n,1);
dzb=zeros(n,1);
%qsback: qs at the upstream node;
%qsit: qs at the node;
%qsfrnt: qs at the downstream node;
%qsdif: spatial difference of qs;
%fI: interfacial exchange fractions;
%dLa: change of activer layer thickness;
%dFi: change of surface fraction;
%delta: change of substrate elevation;
qsback=zeros(n,ng);
qsit=zeros(n,ng);
qsfrnt=zeros(n,ng);
qsdif=zeros(n,ng);
fI=zeros(n,ng);
dLa=zeros(n,1);
dFi=zeros(n,ng);
delta=zeros(n,1);
%zbt: bed elevation at different time;
%st: bed slope at different time;
%qsTt: sediment transport rate per unit width at different time;
%ht: water depth at different time;
zbt=zeros(n,Nt);
st=zeros(n,Nt);
qsTt=zeros(n, Nt);
ht=zeros(n,Nt);
it=1;

%%%%%%Time Processing%%%%%%

while t<Duration
    %%%%%%Flow Hydraulics：Backwater Equation%%%%%%
    %tz=min(t,12000);
    %z(n)=hini+blr*tz/1000;
    z(n)=hini+blr*t/1000;
    h(n)=z(n)-zb(n);
    for ih=1:n-1
        hp=h(n-ih+1)-((zb(n-ih)-zb(n-ih+1))/dx-Cf*qw^2/g/h(n-ih+1)^3)/(1-qw^2/g/h(n-ih+1)^3)*dx;
        h(n-ih)=h(n-ih+1)-0.5*(((zb(n-ih)-zb(n-ih+1))/dx-Cf*qw^2/g/h(n-ih+1)^3)/(1-qw^2/g/h(n-ih+1)^3)+((zb(n-ih)-zb(n-ih+1))/dx-Cf*qw^2/g/hp^3)/(1-qw^2/g/hp^3))*dx;
    end    
    u=qw./h;
    Fr=u./(g.*h).^0.5;
    
    %%%%%%Bed Evolution: Exner equation with active layer formulation%%%%%%
    %%%sediment transport: Sorting Engelund Hansen (An et al., 2021)%%%
    taub=rho.*Cf.*u.^2;
    ustar=u.*(Cf)^0.5;
    shi=taub*ones(1,ng)./rho./R./g./(ones(n,1)*D);
    NA=a.*((ones(n,1)*D)./(Dsg*ones(1,ng))).^b;
    NB=c.*((ones(n,1)*D)./(Dsg*ones(1,ng))).^d;
    Nstar=NA.*(shi).^NB;
    qs=Nstar.*Fi.*((ustar.^3)*ones(1,ng))./R./g./Cf;
    qs=qs*afSEH;
    qsT=sum(qs,2);
    Pstr=qs./(qsT*ones(1,ng));
    Dg_load=2.^(Pstr*Phi')/1000;
    
    %%%evolution of bed elevation%%%
    qsTback=[qsfT;qsT(1:n-1)];
    qsTit=qsT;
    qsTfrnt=[qsT(2:n);2*qsT(n)-qsT(n-1)];
    qsTdif=au*(qsTit-qsTback)+(1-au)*(qsTfrnt-qsTit);
    dzb=-qsTdif/dx*dt*timeyear/(1-p);
    zb=zb+dzb;

    %%%evolution of surface texture%%%
    delta=dzb-dLa;
    fI=((delta<=0)*ones(1,ng)).*Pup+((delta>0)*ones(1,ng)).*(altr*Fi+(1-altr)*Pstr);
    qsback=[qsf;qs(1:n-1,:)];
    qsit=qs;
    qsfrnt=[qs(2:n,:);2*qs(n,:)-qs(n-1,:)];
    qsdif=au*(qsit-qsback)+(1-au)*(qsfrnt-qsit);
    dFi=((-qsdif+qsTdif*ones(1,ng).*fI)/dx*dt*timeyear/(1-p)-(Fi-fI).*(dLa*ones(1,ng)))./(La*ones(1,ng));
    Fi=Fi+dFi;
    Fi=(Fi>0).*Fi;
    Fi=Fi./(sum(Fi,2)*ones(1,ng));

    %%%Stratigraphy storage%%%
    indn=find((delta<=Ls-Lup)&(delta>=-Lup));
    indinc=find(delta>Ls-Lup);
    inddec=find(delta<-Lup);
    Pup(indn,:)=(Pup(indn,:).*(Lup(indn)*ones(1,ng))+fI(indn,:).*(delta(indn)*ones(1,ng)))./((Lup(indn)+delta(indn))*ones(1,ng));
    %Pup=(Pup>0).*Pup;
    %Pup=Pup./(sum(Pup,2)*ones(1,ng));
    Lup(indn)=Lup(indn)+delta(indn);
    if size(indinc,1)>0
        for is=1:size(indinc,1)
            ii=indinc(is);
            Pup(ii,:)=(Pup(ii,:).*Lup(ii)+fI(ii,:).*(Ls-Lup(ii)))./Ls;
            %Pup(ii,:)=(Pup(ii,:)>0).*Pup(ii,:);
            %Pup(ii,:)=Pup(ii,:)./sum(Pup(ii,:),2);
            Store{ii,1}(indup(ii),:)=Pup(ii,:);
            inc=ceil((delta(ii)-(Ls-Lup(ii)))/Ls);
            Store{ii,1}(indup(ii)+1:indup(ii)+inc,:)=ones(inc,1)*fI(ii,:);
            indup(ii)=indup(ii)+inc;
            Pup(ii,:)=fI(ii,:);
            Lup(ii)=delta(ii)-(Ls-Lup(ii))-(inc-1)*Ls;
        end
    end
    if size(inddec,1)>0
        for is=1:size(inddec,1)
            id=inddec(is);
            dec=ceil((-Lup(id)-delta(id))/Ls);
            Store{id,1}(indup(id)-dec+1:indup(id),:)=zeros(dec,ng);
            indup(id)=indup(id)-dec;
            Pup(id,:)=Store{id,1}(indup(id),:);
            Lup(id)=dec*Ls+Lup(id)+delta(id);
        end
    end

    %%%%%%Parameter Update%%%%%%
    t=t+dt;
    i=i+1;
    s=[(zb(1)-zb(2))/dx;(zb(1:n-2)-zb(3:n))/2/dx;(zb(n-1)-zb(n))/dx];
    z=h+zb;
    Dsg=2.^(Fi*Phi')/1000;
    %%%record information evry several years
    if (mod(i,Duration/dt/Nt)==0)
        zbt(:,it)=zb;
        st(:,it)=s;
        qsTt(:,it)=qsT;
        ht(:,it)=h;
        it=it+1;
    end
    
    %figure(1);plot(x,zb,'*r')
end

tend=toc(tstart);
tc=tend/60