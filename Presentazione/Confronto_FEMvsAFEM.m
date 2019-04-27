%% Risoluzione adattiva di un problema del tipo
%
%   -mu''+bu'=f su [a,b]
%   u(a)=u0, u(b)=u1
%
% dove i dati del problema sono contenuti nella struct probdata. La prima
% approssimazione è cercata nello spazio space. La struct method fornisce i
% parametri per l'algoritmo adattivo, come dettagliato nelle functions
% corrispondenti.
clc; clear; close all;
%% Dati del problema
% I dati riportati sono relativi ad un problema di trasporto-diffusione, ma
% si possono modificare a piacere.

probdata.Omega=[0,1];   % Dominio
probdata.m=1;   % Parametri
probdata.b=100;
b=probdata.b; m=probdata.m;
probdata.u0=0;  % Dati al bordo
probdata.u1=1;
probdata.f=@(t) 0.*t; % Termine noto

%% Prima approssimazione
% Calcoliamo una prima approssimazione in uno spazio piccolo per avere
% un'idea della forma della soluzioene e degli errori (stimati) commessi.

% Spazio di approssimazione
space.dim=30;   % con quanti gradi di libertà iniziamo
space.T=linspace(probdata.Omega(1),probdata.Omega(2),space.dim);

% Prima approssimazione
[uh]=solFEM_lin(probdata, space);

% Grafico: produco una figura con la prima approssimazione in modo che
% l'utente possa capire grossomodo con che problema ha a che fare.
figure
plot(space.T,uh)
titolo=sprintf("Approssimazione con %i DoF",space.dim);
title(titolo,'interpreter','latex')

% Stima residuale : calcolo il residuo locale e informo l'utente sul 
% residuo locale e globale nella prima apporssimazione

etaR=LocRes(uh,probdata,space);
etaLoc=max(etaR);
eta=norm(etaR,2);

fprintf("La soluzione determinata ha residuo locale fra %f e %f\n",...
    min(etaR),etaLoc)
fprintf("L'errore globale stimato è %f\n",eta)

stop=input("Si vogliono modificare i parametri della risoluzione adattiva? s/n\n",'s');
if stop=='s'
% Se l'utente vuole dare nuovi parametri alla luce del residuo locale
% appena mostratogli, interrompiamo l'esecuzione per  consentire di
% modificare il codice. 
    return
end

%% Metodo di approssimazione
% In questa sezione vengono raccolti i parametri per il metodo adattivo
% quali il tipo di marker e le soglie. 

method.maxResLoc=Inf; % Massimo errore su ogni elemento della partizione

method.maxRes=1e-5; % Massimo errore globale

method.maxIter=7; % Massimo numero di iterazioni dell'algoritmo adattivo

method.marker='Dor';    % Tipo di marcatore
% method.marker='Max';
method.theta=0.5;   % Parametro del marcatore

method.PreMark=true;    % Sceglie se si effettua premarking
% method.PreMark=false;
method.PreMarkPerc=5;   % Percentuale di elementi da premarcare

%% Risoluzione adattiva
% Vedere la function AFEM.m 
[uh, spacef, err]=AFEM(probdata,space,method);

%% Grafico soluzione finale
subplot(1,2,1)
pbaspect([1,1,1])
hold on
plot(spacef.T,uh) 
plot(spacef.T,min(uh)*ones(size(spacef.T)),'rs')
titolo=sprintf("Approssimazione AFEM con %i DoF",spacef.dim);
title(titolo,'interpreter','latex')

% Informa sull'errore stimato nella soluzioe finale
fprintf("Soluzione trovata con %i gradi di libertà in %i iterazioni\n",...
    spacef.dim, err.NoIter)
fprintf("Abbiamo stimato l'errore locale fra %f e %f\n",min(err.ResLoc),...
    max(err.ResLoc))
fprintf("L'errore globale stimato è %f\n",err.Res)

%% Soluzione FEM
% Spazio di approssimazione
space.dim=spacef.dim;   % eseguo il FEM con lo stesso numeo di DOF
space.T=linspace(probdata.Omega(1),probdata.Omega(2),space.dim);
% Risolvo
[uh]=solFEM_lin(probdata, space);

% Rappresento
subplot(1,2,2)
pbaspect([1,1,1])
hold on
plot(space.T,uh)
plot(space.T,min(uh)*ones(size(space.T)),'rs')
