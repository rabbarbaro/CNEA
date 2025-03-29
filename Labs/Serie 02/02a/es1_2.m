clear
clc
close all

% per risolvere il sistema lineare A*x = b, per quanto sia matematicamente
% corretto, NON facciamo x = inv(A) * b, è estremamente inefficiente:
% utilizziamo invece il comando \, ovvero x = A \ b

A = [50 1 3
    1 6 0
    3 0 1];
n = size(A, 1);

% do un nome alla matrice identità
I = eye(n);

% è possibile determinare i vettori colonna di A^−1 risolvendo n sistemi
% lineari del tipo Av_i = e_i dove e_i sono i vettori della base canonica
% di R^n, con tutti elementi nulli tranne l’elemento i-esimo che è 1

% faccio preventivamente la fattorizzazione LU di A: la matrice è sempre la
% stessa, posso farla una volta sola
[LA, UA, PA] = lu(A);
% prealloco Ainv
Ainv = zeros(n);

for ii = 1:n
    % e_i è uguale alla i-esima colonna della matrice identità
    e = I(:, ii);
    % per risolvere ogni sistema lineare Ax_i = e_i applico fwsub e bksub
    y = fwsub(LA, PA * e);
    x = bksub(UA, y);
    % se non me lo avesse imposto il testo avrei potuto usare x = A / e
    Ainv(:, ii) = x;
    % oppure non lo inizializzo e poi concateno con Ainv = [Ainv x]
    % ma non la posso preallocare, meno efficiente
end

% verifico che A * A^-1 sia la matice identità
A * Ainv;
% calcolo quando mi discosto dalla matrice identità vera
norm(A*Ainv-eye(n));

% verifico la differenza con la matrice inversa trovata con il comando
% inv() di matlab
Ainv - inv(A);