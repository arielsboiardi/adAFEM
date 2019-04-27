function [marked, etaR_rest]=ThreshMark(etaR, threshold)
% Marca tutti gli elementi che hanno residuo locale maggiore della soglia
% threshold.
% Input:    etaR : Vettore dei residui locali
%           threshold : soglia
% Outpu:    marked : struct con 
%               .id : Vettore di indici degli elementi marcati
%               .numel :  numero di elementi marcati
%           etaR_rest : vettore dei residui con gli elementi marcati
%           azzerati

indici=etaR>threshold;
marked.id=find(indici);
marked.numel=nnz(indici);

etaR_rest=etaR;
etaR_rest(marked.id)=0;
end
