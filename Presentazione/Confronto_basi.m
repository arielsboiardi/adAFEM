% %% Confronta la base adattiva con quella uniforme. 
% % Consideriamo il problema 
% %   -mu''+bu'=f su [a,b]
% %   u(a)=u0, u(b)=u1
% %
% % I dati riportati sono relativi ad un problema di trasporto-diffusione.
% 
% clc; clear; close all;
% addpath '..\'
% 
% t=0:0.001:1;
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
% syms z;
% usol=(exp((b/m)*z)-1)/(exp(b/m)-1);
% probdata.uex=matlabFunction(usol);
% 
% %% Metodo di approssimazione
% 
% method.maxResLoc=Inf;
% method.maxRes=0;
% method.maxIter=30;
% method.maxDoF=1000;
% 
% method.marker='Dor';
% % method.marker='Max';
% 
% method.theta=0.5;
% 
% method.PreMark=true;
% % method.PreMark=false;
% method.PreMarkPerc=5;
% 
% %% Inserimento dati dello spazio Xh1
% space.dim=3;
% space.T=linspace(probdata.Omega(1),probdata.Omega(2),space.dim);
% 
% % Visualizzo la base
% subplot(1,2,1)
% plotBase(space)

%% Prima risoluzione === INIZIO ===

[uh, ~]=solFEM_lin(probdata, space);
subplot(1,2,2)
plot(space.T,uh,'LineWidth',2)
hold on
pbaspect([2.5,1,1])
set(gca,'FontName','Calibri Light','FontSize',12);
plot(t,probdata.uex(t),':','LineWidth',2)

%% Stima residuale
% Bisogna a questo punto fornire uno sitmatore residuale nella forma diuna
% function, in quanto lo stimatore è diverso per ogni problema.

etaR=LocRes(uh,probdata,space);

etaLoc=max(etaR);
eta=norm(etaR,2);


%% Marcatura
[marked, etaR]=ThreshMark(etaR, method.maxResLoc);

if method.PreMark
    [marked1,etaR]=PreMark(etaR,method.PreMarkPerc,space);
    marked.id=union(marked.id,marked1.id);
    marked.numel=marked.numel+marked1.numel;
end

switch method.marker
    case 'Dor'
        marked1=DorflerMark(etaR,method.theta,space);
    case 'Max'
        marked1=MaxMark(etaR,method.theta,space);
end

marked.id=union(marked.id,marked1.id);
marked.numel=marked.numel+marked1.numel;

subplot(1,2,1)
plotMarked(marked,space)
subplot(1,2,2)
plotMarked(marked,space)
    
%% Raffinamento dello spazio
space=DyadRef(marked,space);

figure
subplot(1,2,1)
plotBase(space)
close all
space