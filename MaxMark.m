function marked=MaxMark(etaR,theta,space)
% Algoritmo dorfler per marcatura: certa riduzione dell'errore.
% Input:    etaR: vettore con i residui locali 
%           theta: soglia 0<theta<1
%           space: struct con le informazioni dello spazio (vedi AFEM.m)
% Output:   marked: struct con i seguenti campi
%              - id: indici degli elementi marcati
%              - numel: numero di elementi marcati
N=space.dim-1;
idMarked=[];
NoMarked=0;
eta=max(etaR); 
for idK=1:N
    if etaR(idK)>=theta*eta
        idMarked=[idMarked,idK];
        NoMarked=NoMarked+1;
    end
end
marked.id=idMarked;
marked.numel=NoMarked;
end