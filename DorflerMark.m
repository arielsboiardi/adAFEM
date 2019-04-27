function marked=DorflerMark(etaR,theta,space)
% Algoritmo di marcatura di Dorfler
% Input:    etaR: vettore con i residui locali 
%           theta: soglia 0<theta<1
%           space: struct con le informazioni dello spazio di splines
% Output:   marked: struct con i seguenti campi
%              - id: indici degli elementi marcati
%              - numel: numero di elementi marcati
N=space.dim-1;
idMarked=[];
NoMarked=0;
Sigma=0; % Rappresenta (errore dovuto agli elementi selezionati)^2
err_glob=sum(etaR.^2); % Rappresenta (errore globale)^2

while Sigma<theta*err_glob
    Tt=setdiff([1:N],idMarked);
    eta=max(etaR(Tt));
    for K=Tt
        if etaR(K)==eta
            idMarked=[idMarked,K];
            NoMarked=NoMarked+1;
            Sigma=Sigma+etaR(K)^2;
        end
    end
end

marked.id=idMarked;
marked.numel=NoMarked;
end