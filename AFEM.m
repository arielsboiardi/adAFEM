function [uh, spacef, err]=AFEM(probdata, space, method)
% Realizza l'algoritmo degli elementi finiti adattivo per un problema della
% forma 
%
%   -m u'' + b u' = f su [a,b]
%   u(a)=u0, u(b)=uN
%
% dove i dati del problema sono contenuti nella struct probdata. La prima
% approssimazione è cercata nello spazio space. La struct method fornisce i
% parametri per l'algoritmo adattivo. 
% 
% Input:    probdata : struct che rappresenta i dati del problema
%               .Omega : Vettore 1x2 con estremi a,b dell'intervallo
%               .m, .b : Parametri dell'equazione
%               .u0, .uN : Valori al bordo
%               .f : Termine noto dell'equazione differenziale
%
%           space : struct contenente le informazioni dello spazio di
%                   approssimazione iniziale
%               .dim : Dimensione dello spazio di approssimazione
%               .T : Vettore dei nodi della partizione
%           method : struct contenente i parametri dell'algoritmo adattivo
%               .maxIter : Massimo numero di iterazioni
%               .maxDoF : Massimo numero di gradi di libertà
%               .maxRes : Massimo errore stimato mediante residuo accettato
%               .maxResLoc : Massimo residuo accettato su ogni elemento
%                            della partizione
%               .marker : Stringa == 'Dor' per il marcatore di Dorfler
%                                 == 'Max' per il criterio del massimo
%               .theta : Parametro per il marcatore
%               .PreMark : Booleano che indica se il premarking sia attivo
%               .PreMarkPerc : Percentuale di elementi da premarcare
%
% Output:   uh : Coordinate della soluzione approssimata rispetto alla base
%                dello spazio descritto da spacef
%           spacef : struct uguale a space contenente le informazioni dello
%                    spazio di approssimazione costruito con il metodo 
%                    adattivo
%           err : struct con le informazioni relative all'errore stimato
%               .Res : Residuo globale dell'ultima approssimazione
%               .ResLoc : Residuo locale dell'ultima approssimazione
%               .NoIter : Numero di iterazioni richiesto

% Controllo che la partizione su cui si costruisce lo spazio space sia
% compatibile con il problema assegnato.
if (space.T(1)~=probdata.Omega(1)) || (space.T(end)~=probdata.Omega(2))
    error("La partizione space.T non è valida per il problema assegnato")
end

% Inizalizzazioni
eta=Inf;    % Per essere certi che almeno venga svolta un'iterazione
space_ref=space; % Lo spazio inizia come quello dato, poi viene modificato

for kdx=0:method.maxIter
    % Aggiornamento dati 
    spacef=space_ref;
    
    % Approssimazione
    uh=solFEM_lin(probdata,spacef); 
    
    % Stima il residuo locale su tutti gli elementi della partizione
    etaR=LocRes(uh,probdata,spacef);
    eta=norm(etaR,2);   % Residuo globale
    
    err.ResLoc=etaR;
    err.Res=eta;
    
    % Test per uscita a soglia del residuo globale
    if eta<=method.maxRes
        break
    end
    
    % Marcatore a soglia 
    % Marca per il raffinamento tutti gli elementi con errore stimato
    % maggiore di sigma=method.maxResLoc
    [marked, etaR]=ThreshMark(etaR, method.maxResLoc);
    
    % Premark
    % marca il (method.PreMarkPerc)% degli elementi con errore massimo se
    % l'utente lo ha chiesto
    if method.PreMark
        [marked1,etaR]=PreMark(etaR,method.PreMarkPerc,spacef);
        % Unisco gli elementi fino ad ora marcati
        marked.id=union(marked.id,marked1.id);   
        marked.numel=marked.numel+marked1.numel;
    end
    
    % Marcatori
    switch method.marker
        case 'Dor'
            marked1=DorflerMark(etaR,method.theta,spacef);
        case 'Max'
            marked1=MaxMark(etaR,method.theta,spacef);
    end
    
    % Unisco tutti gli elementi marcati 
    marked.id=union(marked.id,marked1.id);
    marked.numel=marked.numel+marked1.numel;
    
    % Raffinamento dello spazio di approssimazione
    space_ref=DyadRef(marked,spacef);
end
% Nel caso di uscita a soglia restituisco il numero di iterazioni che è
% stato necessario.
err.NoIter=kdx; 
end
    
    