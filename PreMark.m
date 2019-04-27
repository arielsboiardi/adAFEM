function [marked,etaR_rest]=PreMark(etaR,eps,space)
% Premarker all'eps%
% Input:    etaR: vettore col residuo locale
%           eps: percentuale degli elementi da premarcare
% Output:   marked: struct con
%              - id: indici degli elementi marcati
%              - numel: numero di elementi marcati
%           etaR_rest: vettore dei residui azzerato sugli elementi
%           premarcati, in modo che non vengano marcati successivamente 

N=space.dim-1; % numero elementi di etaR
etaR_rest=etaR; % lo salvo

NoMarked=fix(N*eps*0.01); % NoMarked è il p% degli elementi (intero)
if nnz(etaR)>NoMarked
    [~,idMarked]=sort(etaR,'descend');
    idMarked=idMarked(1:NoMarked);
    etaR_rest(idMarked)=0; % Azzero il residuo corrispondente agli elementi 
                           % marcati
    marked.id=idMarked;
    marked.numel=NoMarked;
else
    marked.id=[];
    marked.numel=0;
end
end
    