%% Convergenza in norma L2 del metodo FEM e AFEM

clc; clear; close all;

%% Costruiamo i dati del problema data la soluzione
probdata.Omega=[0,1];  % Dominio
probdata.m=1;   % Parametri
probdata.b=500;
b=probdata.b; m=probdata.m;
probdata.u0=0;  % Dati al bordo
probdata.u1=1;

% k=15;
% syms z 
% u=(cos(k*z)).^1;
% du=diff(u,z);
% ddu=diff(du,z);
% f=-m*ddu+b*du;

probdata.f=@(t) 0.*t;

probdata.uex=@(x) (exp((b/m)*x)-1)/(exp(b/m)-1);
% probdat.duex=matlabFunction(du);

%% Inserimento dati dello spazio Xh1
space.dim=100;
space.T=linspace(probdata.Omega(1),probdata.Omega(2),space.dim);

%% Prima risoluzione
[uh, Uh]=solFEM_lin(probdata, space);

% Grafico
plot(space.T,uh)
titolo=sprintf("Approssimazione con %i DoF",space.dim);
title(titolo,'interpreter','latex')

%% Stima residuale
% Bisogna a questo punto fornire uno sitmatore residuale nella forma diuna
% function, in quanto lo stimatore è diverso per ogni problema.

etaR=LocRes(uh,probdata,space);

etaLoc=max(etaR);
eta=norm(etaR,2);

fprintf("La soluzione determinata ha residuo locale fra %f e %f\n",...
    min(etaR),etaLoc)
fprintf("L'errore globale stimato è %f\n",eta)

stop=input("Si vogliono cambiare i parametri? s/n\n",'s');
if stop=='s'
    return
end

%% Errore L2
kdx=1;

DOF(kdx)=space.dim;

errsq=@(t) (Uh(t)-probdata.uex(t)).^2;
Q=integral(errsq,0,1);
ERRL2(kdx)=sqrt(Q);
ETA(kdx)=eta;

%% Metodo di approssimazione

method.maxResLoc=Inf;
method.maxRes=0;
method.maxIter=1000;
method.maxDoF=1000;

method.marker='Dor';
% method.marker='Max';

method.theta=0.5;

method.PreMark=false;
% method.PreMark=false;
method.PreMarkPerc=5;

%% Ciclo adattivo
while (eta>method.maxRes) && (space.dim<method.maxDoF)
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
    
    %% Raffinamento dello spazio
    space=DyadRef(marked,space);
    
    %% Risoluzione
    [uh, Uh]=solFEM_lin(probdata, space);
    
    %% Stima
    etaR=LocRes(uh,probdata,space);
    eta=norm(etaR,2);
    etaLoc=max(etaR);
    
    %% Errore L2
    DOF(kdx)=space.dim;
    
    errsq=@(t) (Uh(t)-probdata.uex(t)).^2;
    Q=integral(errsq,0,1);
    ERRL2(kdx)=sqrt(Q);
    ETA(kdx)=eta;
    kdx=kdx+1;
end

%% Soluzione finale
t=probdata.Omega(1):0.1/space.dim:probdata.Omega(2);
plot(space.T,uh)
hold on
plot(space.T,0*ones(size(space.T)),'rs')
titolo=sprintf("Approssimazione AFEM con %i DoF",space.dim);
title(titolo,'interpreter','latex')

fprintf("Soluzione trovata con %i gradi di libertà\n",space.dim)
fprintf("Abbiamo stimato l'errore locale fra %f e %f\n",min(etaR),etaLoc)
fprintf("L'errore globale stimato è %f\n",eta)

%% Confronto col FEM

for kdx=1:numel(DOF)
    DoF=DOF(kdx);
    space.dim=DoF;
    space.T=linspace(probdata.Omega(1),probdata.Omega(2),space.dim);
    [uh, Uh]=solFEM_lin(probdata, space);

    %% Stima errore
    etaR=LocRes(uh,probdata,space);
    eta=norm(etaR,2);
    ETA_FEM(kdx)=eta;

    %% Errore L2 FEM
    errsq=@(t) (Uh(t)-probdata.uex(t)).^2;
    Q=integral(errsq,0,1);
    ERRL2_FEM(kdx)=sqrt(Q);
    ETA(kdx)=eta;
end

% Rappresentazione 
t=probdata.Omega(1):0.1/space.dim:probdata.Omega(2);
figure
plot(space.T,uh)
hold on
plot(space.T,0*ones(size(space.T)),'rs')
titolo=sprintf("Approssimazione FEM con %i DoF",space.dim);
title(titolo,'interpreter','latex')

%% Convergenza
figure
hold on
set(gca,'Yscale','log')
set(gca,'Xscale','log')
plot(DOF,ERRL2_FEM,'kd-')
plot(DOF,ERRL2,'ks-')
legend("FEM","AFEM")
% titolo=strcat("Errore $L^2$ ($k=",num2str(k),"$)");
% title(titolo, 'interpreter','latex')
