function D = qrbasic_shift(A,tol,nmax)
%
% D = qrbasic_shift(A,tol,nmax)
%
% Metodo delle iterazioni QR per l'approssimazione degli n autovalori di 
% una matrice A quadrata a valori reali.
%
% Parametri di ingresso:
% A       Matrice a valori reali quadrata (n x n)
% tol     Tolleranza sul criterio d'arresto
% nmax    Numero massimo di iterazioni
%
% Parametri di uscita:
% D       Vettore colonna (n x 1) contenente gli n autovalori approssimati
%
%                                         Politecnico di Milano, 04/04/2024
%

[n,m]=size(A);

mu = 0;
niter = 0;
test = max(max(abs(tril(A,-1))));
while niter < nmax && test > tol
  [Q,R]=qr(A - mu*eye(n));    
  A = R*Q + mu*eye(n);
  mu = A(n, n);
  niter = niter + 1;
  test = max(max(abs(tril(A,-1))));
end
if niter > nmax
 fprintf(['Il metodo non converge nel massimo',...
             ' numero di iterazioni permesso']);
else
 fprintf('Il metodo converge in %d iterazioni\n',niter)
end
D = diag(A);
return
