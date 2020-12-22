clear;
clc;
tstart=tic;

%%%%%%%%%%%%Genaralized Simulation%%%%%%%%%%%%

Ltilde=10;
dxtilde=0.1;
n=Ltilde/dxtilde+1;
Duration=30;
dttilde=1e-3;
Frn=0.8;
Bl=1;
au=1;
Nt=15;


%%%%%%Initial condition%%%%%%
%xtilde: dimensionless distance;
%Htilde: dimensionless water depth;
%etadtilde: dimensionless deviatoric bed elevation;
%Sdtilde: dimensionless deviatoric bed slope;
%qstilde: dimensionless sediment transport rate per unit width;
%ztilde: dimensionless water surface at junction;
%ztilde_ini: initial ztilde;
xtilde=linspace(0,Ltilde,n)';
Htilde=ones(n,1);
etadtilde=zeros(n,1);
Sdtilde=[(etadtilde(1)-etadtilde(2))/dxtilde;(etadtilde(1:n-2)-etadtilde(3:n))/2/dxtilde;(etadtilde(n-1)-etadtilde(n))/dxtilde];
qstilde=zeros(n,1);
ztilde_ini=etadtilde(n)+Htilde(n);
ttilde=0;
i=0;
%ttilde_incip: incipient ttilde for blocked valley lake;
%xtilde_incip: incipient xtilde for blocked valley lake;
ttilde_incip=0;
xtilde_incip=0;
%qstildeback: qstilde at the upstream node;
%qstildeit: qstilde at the node;
%qstildefrnt: qstilde at the downstream node;
%qstildedif: spatial difference of qstilde;
%detadtilde: change of etadtilde;
qstildeback=zeros(n,1);
qstildeit=zeros(n,1);
qstildefrnt=zeros(n,1);
qstildedif=zeros(n,1);
detadtilde=zeros(n,1);
%etadtilde_store: etadtilde at different time;
%Sdtilde_store: Sdtilde at different time;
%qstilde_store: qstilde at different time;
%Htilde_store: Htilde at different time;
etadtilde_store=zeros(n,Nt);
Sdtilde_store=zeros(n,Nt);
qstilde_store=zeros(n, Nt);
Htilde_store=zeros(n,Nt);
it=1;

%%%%%%Time Processing%%%%%%

while ttilde<Duration
    %%%Flow Hydraulics£ºBackwater Equation%%%
    ztilde=ztilde_ini+Bl*ttilde;
    Htilde(n)=ztilde-etadtilde(n);
    for ih=1:n-1
        Htilde_p=Htilde(n-ih+1)-(1+(etadtilde(n-ih)-etadtilde(n-ih+1))/dxtilde-1/Htilde(n-ih+1)^3)/(1-Frn^2/Htilde(n-ih+1)^3)*dxtilde;
        Htilde(n-ih)=Htilde(n-ih+1)-0.5*((1+(etadtilde(n-ih)-etadtilde(n-ih+1))/dxtilde-1/Htilde(n-ih+1)^3)/(1-Frn^2/Htilde(n-ih+1)^3)+(1+(etadtilde(n-ih)-etadtilde(n-ih+1))/dxtilde-1/Htilde_p^3)/(1-Frn^2/Htilde_p^3))*dxtilde;
    end
    
    %%%sediment transport: modified Engelund Hansen%%%
    qstilde=Htilde.^(-5);
    
    %%%Exner equation%%%
    qstildeback=[1;qstilde(1:n-1)];
    qstildeit=qstilde;
    qstildefrnt=[qstilde(2:n);2*qstilde(n)-qstilde(n-1)];
    qstildedif=au*(qstildeit-qstildeback)+(1-au)*(qstildefrnt-qstildeit);
    detadtilde=-qstildedif/dxtilde*dttilde;
    etadtilde=etadtilde+detadtilde;
    
    %%%parameter update%%%
    ttilde=ttilde+dttilde;
    i=i+1;
    Sdtilde=[(etadtilde(1)-etadtilde(2))/dxtilde;(etadtilde(1:n-2)-etadtilde(3:n))/2/dxtilde;(etadtilde(n-1)-etadtilde(n))/dxtilde];
    
    %%%incipient of blocked valley lake%%%
    if ttilde_incip==0
       Sdtildedif=[0;Sdtilde(2:n)-Sdtilde(1:n-1)];
       ind=find(Sdtildedif>0.001);
       ttilde_incip=(1-isempty(ind))*ttilde;
       xtilde_incip=find(Sdtildedif>0.001,1)*dxtilde-dxtilde;
    end
    
    %%%store data%%%
    if (mod(i,Duration/dttilde/Nt)==0)
        etadtilde_store(:,it)=etadtilde;
        Sdtilde_store(:,it)=Sdtilde;
        qstilde_store(:,it)=qstilde;
        Htilde_store(:,it)=Htilde;
        it=it+1;
    end
    
end

tend=toc(tstart);
tc=tend/60