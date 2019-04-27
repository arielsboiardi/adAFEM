function refined_space=DyadRef(marked,space)
% Raffinamento diadico degli elementi della mesh dello spazio space marcati
% in marked.
% Input:    marked: struct (si veda una function di marcatura)
%           space: informazioni dello spazio (si veda AFEM.m)
% Output:   refined_space: struct basata su space che rappresenta lo spazio
%           raffinato.

refined_space=space; % iniziamo dallo spazio dato
T=space.T; % la partizione su cui è basato è 

% Siccome costruisco una riga controllo che la partizione sia in formato
% riga. Se non la è la metto in riga.
temp=size(T);
if temp(1)<temp(2)
    T=reshape(T,1,[]);
end
clear temp

for jdx=1:marked.numel
    k=marked.id(jdx);
    h=(T(k+1)-T(k));
    if h>0
        T=[T(1:k),T(k)+0.5*h,T(k+1:end)];
        refined_space.dim=refined_space.dim+1;
        marked.id=marked.id+1;
    end
end
refined_space.T=T;
refined_space.dim=numel(T);
end