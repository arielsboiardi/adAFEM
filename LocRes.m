function etaR=LocRes(uh,probdata,space)
% Stimatore residuale per problemi della forma
%   -mu''+bu'=f su [a,b] &   u(a)=u0, u(b)=u1
%
% Input:    uh : Valori nodali dell'approssimazione
%           probdata : struct con i dati del problema (vedi AFEM.m)
%           space : struct con i dati dello spazio (vedi AFEM)
% Output:   etaR : vettore con i residui locali

b=probdata.b;   % Leggo i dati
m=probdata.m;
f=probdata.f;

T=space.T;      % Imposto lo spazio
N=space.dim-1;
H=diff(T);

% Se il termine noto è nullo abbiamo un'espressione chiusa. Altrimenti
% busogna integrare...
t=probdata.Omega(1):0.0001:probdata.Omega(2);
nnull=nnz(f(t));

etaR=zeros(1,N);
for kdx=1:N
    if nnull==0
        etaR(kdx)=b^2*sqrt(H(kdx))*abs(uh(kdx+1)-uh(kdx));
    else
        fun=@(t) (f(t)-b^2*(1/H(kdx))*(uh(kdx+1)-uh(kdx))).^2;
        Q=integral(fun,T(kdx),T(kdx+1));
        etaR(kdx)=H(kdx)*sqrt(Q);
    end
end
end

