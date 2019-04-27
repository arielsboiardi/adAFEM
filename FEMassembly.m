function A=FEMassembly(probdata,space)
% Costruisce la matrice di rigidezza del problema probdata nello spazio
% space.

m=probdata.m;
b=probdata.b;
T=space.T;
N=space.dim-1;
H=diff(T);  % Passi della partizione

A=zeros(N-1,N-1);
for idx=1:N-1
    if idx~=1
        A(idx,idx-1)=-m/H(idx)-b/2;
    end
    A(idx,idx)=m*(1/H(idx)+1/H(idx+1));
    if idx~=N-1
        A(idx,idx+1)=-m/H(idx+1)+b/2;
    end
end
end