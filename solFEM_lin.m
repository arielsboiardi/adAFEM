function [uh,Uh]=solFEM_lin(probdata,space)
% Determina la soluzione approssimata nello spazio X^1 definito da space
% con dati probdata del problema 
%
%   -m u'' + b u' = f su [a,b]
%   u(a)=u0, u(b)=uN
%
% Input:    probdata : struct che rappresenta i dati del problema
%               .Omega : Vettore 1x2 con estremi a,b dell'intervallo
%               .m, .b : Parametri dell'equazione
%               .u0, .uN : Valori al bordo
%               .f : termine noto dell'equazione differenziale
%
%           space : struct contenente le informazioni dello spazio di
%                   approssimazione
%               .dim : Dimensione dello spazio di approssimazione
%               .T : Vettore dei nodi della partizione
%
% Output:   uh: vettore con i valori nodali di uh
%           Uh: soluzione Uh come function handle (facoltativo)

% Calcolo la soluzione
A=FEMassembly(probdata,space);
F=FEMload(probdata, space);
uh=A\F;

% completo la soluzione con i dati di bordo
uh=[probdata.u0, uh', probdata.u1];

if nargout==2
% Se l'utente ne ha bisogno fornisco già una versione funzionale della
% soluzione
    Uh=@(s) interp1(space.T,uh,s);
end
end
