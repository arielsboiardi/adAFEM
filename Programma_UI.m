%% Risoluzione adattiva di un problema di trasporto-diffusione
% Questo codice serve per raccogliere un esempio ocmpleto e generale ed
% implicitamente definire le struct che useremo nei passi successivi.
%
% Abbiamo pensato di scrivere tutto usando delle struct invece che lasciare
% la variabili isolate in modo da migliorare la leggibilità dei codici e
% aumentare la flessibilità del lavoro.
clc

%% Inserimento dei dati del problema
fprintf("Vuoi dare dei nuovi dati o usare dati già esistenti?\n")
fprintf("    a. Nuovi dati \n    b. Dati da file\n")
if exist('probdata')
    fprintf("    c. Dati in workspace\n")
end
caricamento=input("Scelta: ",'s');

switch caricamento
    case 'a'
        clear
        fprintf("Risolvo un poblema del tipo:\n -mu''+bu'=f su [a,b]\n u(0)=u0, u(1)=u1\n")
        fprintf("Inserire i dati\n")
        
        probdata.Omega=input("    Estremi dell'intervallo: [a, b]= ");
        
        probdata.m=input("    m= ");
        probdata.b=input("    b= ");
        
        probdata.u0=input("    u(a)= ");
        probdata.u1=input("    u(b)= ");
        
        probdata.f=input("    f= ");
    case 'b'
        datafile=input("Fornire il nome del file contenente i dati");
        load(datafile);
    case 'c'
        if exist('probdata')
            fprintf("D'accordo\n");
        else
            restart=input("Non ci sono dati, ricominciamo? s/n \n",'s');
            switch restart
                case 's'
                    Programma
                case 'n'
                    return
            end
        end
end

chiudi=input("Vuoi cancellare le figure? s/n \n",'s');
switch chiudi
    case 's'
        close all
    case 'n'
        fprintf("Mantengo le figure")
        hold on
end

%% Inserimento dati dello spazio Xh1
fprintf("Calcoliamo una prima soluzione con FEM lineare.\n")

space.dim=input("Con quanti DoF inizio? ");
space.T=linspace(probdata.Omega(1),probdata.Omega(2),space.dim);

%% Prima risoluzione
[uh, Uh]=solFEM_lin(probdata, space);

% Grafico
t=probdata.Omega(1):0.1/space.dim:probdata.Omega(2);
plot(space.T,uh)
titolo=sprintf("Approssimazione con %i DoF",space.dim);
title(titolo,'interpreter','latex')

%% Stima residuale

etaR=LocRes(uh,probdata,space);

etaLoc=max(etaR);
eta=norm(etaR,2);

%% Metodo di approssimazione
if exist('method')
    scelta=input("Vuoi stabilire nuovi parametri per la risoluzione? s/n \n",'s');
else
    scelta='s';
end

switch scelta
    case 's'
        fprintf("La soluzione determinata ha residuo locale fra %f e %f\n", min(etaR),etaLoc)
        fprintf("L'errore globale stimato è %f\n",eta)
        
        method.maxResLoc=input("Quale è il masismo residuo locale accettabile?\n");
        %         fprintf("Accettermo soluzioni con residuo locale minore di %d\n",...
        %             method.maxResLoc)
        method.maxRes=input("Quale è il masismo residuo globale accettabile?\n");
        
        method.maxDoF=input("Fino a quanti gradi di libertà posso arrivare?\n");
        
        fprintf("Che marcatore vuoi usare?\n")
        fprintf("    a. Dorfler/Equilibrio\n    b. Massimo\n")
        scelta=input("Scelta: ",'s');
        
        switch scelta
            case 'a'
                method.marker='Dor';
                fprintf("Utilizziamo il marcatore Dorfler")
            case 'b'
                method.marker='Max';
                fprintf("Utilizziamo il marcatore del massimo")
        end
        method.theta=input(" con soglia theta=");
        
        scelta=input("Vuoi usare un premarker? s/n \n",'s');
        switch scelta
            case 's'
                method.PreMark=true;
                method.PreMarkPerc=input("Che percentuale di elementi premarchiamo?\n");
            case 'n'
                method.PreMark=false;
        end
        
    case 'n'
        fprintf("Continuo a fare come ho fatto fino ad ora")
end

%% Ciclo adattivo
% while (etaLoc>method.maxResLoc) && ...
%         (eta>method.maxRes) && 
while    (space.dim<method.maxDoF)
    %% Marcatura
    [marked, etaR]=ThreshMark(etaR, method.maxResLoc);
    
    if method.PreMark
        [marked1,etaR]=PreMark(etaR,method.PreMarkPerc,space);
        marked.id=union(marked.id,marked1.id);
        marked.numel=marked.numel+marked1.numel;
    end
    
    switch method.marker
        case 'Dorfler'
            marked1=DorflerMark(etaR,method.theta,space);
        case 'massimo'
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
end

%% Soluzione finale
plot(space.T,uh)
hold on

vedinodi=input("vuoi vedere i nodi? s/n \n",'s');
switch vedinodi
    case 's'
        plot(space.T,min(abs(Uh(t)))*ones(size(space.T)),'r.')
        legend("$u_h$","Nodi",'interpreter','latex')
    case 'n'
end
titolo=sprintf("Approssimazione AFEM con %i DoF",space.dim);
title(titolo,'interpreter','latex')

fprintf("Soluzione trovata con %i gradi di libertà\n",space.dim)
fprintf("Abbiamo stimato l'errore locale fra %f e %f\n",min(etaR),etaLoc)
fprintf("L'errore globale stimato è %f\n",eta)


