% %% Epprossimazione adattiva ed errore puntuale
% addpath('..\');
% clc; clear; close all;
% 
% %% Costruiamo i dati del problema data la soluzione
% probdata.Omega=[0,1];  % Dominio
% probdata.m=1;   % Parametri
% probdata.b=100;
% b=probdata.b; m=probdata.m;
% probdata.u0=0;  % Dati al bordo
% probdata.u1=1;
% 
% probdata.f=@(t) 0.*t;
% 
% probdata.uex=@(x) (exp((b/m)*x)-1)/(exp(b/m)-1);
% 
% % Spazio iniziale
% space.dim=30;
% space.T=linspace(probdata.Omega(1),probdata.Omega(2),space.dim);
% 
% %% Metodo adattivo 
% method.maxResLoc=10000;
% method.maxRes=0;
% method.maxIter=3000;
% method.maxDoF=1000;
% 
% method.marker='Dor';
% % method.marker='Max';
% 
% method.theta=0.9;
% 
% method.PreMark=true;
% % method.PreMark=false;
% method.PreMarkPerc=5;

%% Approssimazione === INIZIO ===
[uh, Uh]=solFEM_lin(probdata, space);
% 
% %Grafici
% err=@(t) abs(probdata.uex(t)-Uh(t));
% s=0:0.0001:1;
% 
% subplot(1,2,1)
% plot(space.T,uh,'LineWidth',2)
% pbaspect([1,1,1]);
% hold on
% 
% subplot(1,2,2)
% plot(s,err(s),'LineWidth',2);
% set(gca,'YScale','log')
% hold on

%% Stima residuale
etaR=LocRes(uh,probdata,space);

etaLoc=max(etaR);
eta=norm(etaR,2);

%% Marcatura
[marked, etaR]=PreMark(etaR, method.PreMarkPerc,space);

switch method.marker
    case 'Dor'
        marked1=DorflerMark(etaR,method.theta,space);
    case 'Max'
        marked1=MaxMark(etaR,method.theta,space);
end

marked.id=union(marked.id,marked1.id);
marked.numel=marked.numel+marked1.numel;

%% Raffinamento dello spazio
space=DyadRef(marked,space);
space

