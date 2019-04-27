%% Confronto errore-residuo
% clear all
addpath '..\'

%% Costruiamo i dati del problema data la soluzione
probdata.Omega=[0,1];  % Dominio
probdata.m=1;   % Parametri
probdata.b=10;
b=probdata.b; m=probdata.m;
probdata.u0=0;  % Dati al bordo
probdata.u1=1;

probdata.f=@(t) 0.*t;

syms z;
usol=(exp((b/m)*z)-1)/(exp(b/m)-1);
dusol=diff(usol,z);
probdata.uex=matlabFunction(usol);
probdata.duex=matlabFunction(dusol);

%% Inserimento dati dello spazio Xh1
% space.dim=31;
% space.T=linspace(probdata.Omega(1),probdata.Omega(2),space.dim);

%% Prima risoluzione === INIZIO ===
T=space.T;
H=diff(T);

[uh, Uh]=solFEM_lin(probdata, space);

err=@(t) abs(probdata.uex(t)-Uh(t));
t=0:0.0001:1;
plot(t,err(t))
hold on

%% Stima residuale

etaR=LocRes(uh,probdata,space);

eta=norm(etaR,2);

%% Errore in norma H1
for kdx=1:space.dim-1
    errsq=@(t) (probdata.uex(t)-((uh(kdx)/H(kdx))*(t-T(kdx))...
        +(uh(kdx+1)/H(kdx))*(T(kdx+1)-t))).^2;
    derrsq=@(t) (probdata.duex(t)-1/H(kdx)*(uh(kdx+1)-uh(kdx))).^2;
    Q1=integral(errsq,T(kdx),T(kdx+1));
    Q2=integral(derrsq,T(kdx),T(kdx+1));
    ErrH1(kdx)=sqrt(Q1)+sqrt(Q2);
%     ErrH1(kdx)=sqrt(Q1);
end

plot(H/2+T(1:end-1),ErrH1,'ko')
hold on
plot(H/2+T(1:end-1),etaR,'ks')
set(gca,'YScale','log')

