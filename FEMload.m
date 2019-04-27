function F=FEMload(probdata,space)
% Costruisce il vettore di carico per il probelma definito in probdata
% nello spazio space.

m=probdata.m;   % Estraggo i dati
b=probdata.b;
T=space.T;
N=space.dim-1;

H=diff(T); % Passi della partizione
fn=probdata.f;
F=zeros(space.dim-2, 1);    % Inizalizzo

for kdx=1:N-1
    fn1=@(t) fn(t).*(t-T(kdx)); % Contributi del termine noto
    fn2=@(t) fn(t).*(T(kdx+2)-t);
    Q1=integral(fn1,T(kdx),T(kdx+1));
    Q2=integral(fn2,T(kdx+1),T(kdx+2));
    F(kdx)=Q1/H(kdx)+Q2/H(kdx+1);
end
F(1)=F(1)+(m/H(1) + 0.5*b)*probdata.u0; % Contributo dei dati al bordo
F(N-1)=F(N-1)+(m/H(N-1) - 0.5*b)*probdata.u1;
end