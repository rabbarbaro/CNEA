clear
clc
close all

% definisco la dimensione del sistema e costruisco la matrice e il vettore
% dei termini noti come da testo
n = 4;
A = diag(10.^(0:n-1)) - diag(3*ones(n-1, 1), 1) - diag(2*ones(n-2, 1), 2) ...
    - diag(ones(n-3, 1), 3) - diag(3*ones(n-1, 1), -1) ...
    - diag(2*ones(n-2, 1), -2) - diag(ones(n-3, 1), -3);
b = [-5, 2, 92, 994]';

%% 1

% dalla teoria:
% grad: ||x-x_k||_A <= d^k * ||x-x0||
%       d = (K(A)-1) / (K(A)+1)
% grad conj: ||x-x_k||_A <= (2*c^k)/(1+c^(2*k)) * ||x-x0||_A
%            c = (sqrt(K(A))-1) / (sqrt(K(A))+1)

% dopo 100 iterazioni
k = 100;

% calcolo il numero di condizionamento spettrale di A
K_A = max(eig(A)) / min(eig(A));

% per il gradiente calcolo d^k
d = (K_A - 1) / (K_A + 1);
fatt_abb_erro_grad_100 = d^k;
fprintf("Il fattore di abbattimento dell'errore per metodo del gradiente dopo %d iterazioni: %f\n", ...
    k, fatt_abb_erro_grad_100);

% per il gradiente coniugato calcolo (2*c^k)/(1+c^(2*k))
c = (sqrt(K_A)-1) / (sqrt(K_A)+1);
fatt_abb_erro_grad_conj_100 = (2*c^k) / (1+c^(2*k));
fprintf("Il fattore di abbattimento dell'errore per metodo del gradiente coniugato dopo %d iterazioni: %f\n", ...
    k, fatt_abb_erro_grad_conj_100);

%% 2

% eseguiamo a mano l'algoritmo del gradiente non precondizionato per le
% prime 2 iterazioni (oltre alla 0)

% inizializziamo guess iniziale, il residuo iniziale e alpha al passo 0
x0 = zeros(size(b));
r0 = b - A*x0;
alpha0 = (r0' * r0) / (r0' * A * r0);

% aggiorniamo la soluzione, il residuo e alpha al passo 1
x1 = x0 + alpha0 * r0;
r1 = r0 - alpha0 * A * r0;
alpha1 = (r1' * r1) / (r1' * A * r1);

% aggiorniamo soluzione e residuo al passo 2
x2 = x1 + alpha1 * r1;
r2 = r1 - alpha1 * A * r1;

% nel metodo del gradiente non precondizionato le direzioni di discesa sono
% i residui, devo calolare l'angolo fra residui successivi
% angolo fra v1 e v2: arccos((v(1)/||v(1)||)' * (v(2)/||v(2)||))
angolo1 = acos((r0 / norm(r0))' * (r1 / norm(r1)));
angolo2 = acos((r1 / norm(r1))' * (r2 / norm(r2)));

fprintf("Angolo fra r0 ed r1: %f\n", rad2deg(angolo1));
fprintf("Angolo fra r1 ed r2: %f\n", rad2deg(angolo2));

%% 3

% eseguiamo a mano l'algoritmo del gradiente coniugato (non
% precondizionato) per le prime 2 iterazioni (oltre alla 0)

% inizializziamo guess iniziale, il residuo iniziale e la direzione di
% dicesa al passo 0
x0 = zeros(size(b));
r0 = b - A*x0;
p0 = r0;

% applico l'algoritmo per k = 0 (alpha0, x1, r1, beta0, p1)
alpha0 = (p0' * r0) / (p0' * A * p0);
x1 = x0 + alpha0 * p0;
r1 = r0 - alpha0 * A * p0;
beta0 = (p0' * A * r1) / (p0' * A * p0);
p1 = r1 - beta0 * p0;

% applico l'algoritmo per k = 1 (alpha1, x2, r2, beta1, p2)
alpha1 = (p1' * r1) / (p1' * A * p1);
x2 = x1 + alpha1 * p1;
r2 = r1 - alpha1 * A * p1;
beta1 = (p1' * A * r2) / (p1' * A * p1);
p2 = r2 - beta1 * p1;

% nel metodo del gradiente coniugato non precondizionato le direzioni di
% discesa sono tutte A-ortogonali le une alle altre
% angolo fra v1 e v2: arccos((v(1)/||v(1)||)' * (v(2)/||v(2)||))
angolo1 = acos((p0 / norm(p0))' * (p1 / norm(p1)));
angolo2 = acos((p1 / norm(p1))' * (p2 / norm(p2)));

fprintf("Angolo fra p0 ed p1: %f\n", rad2deg(angolo1));
fprintf("Angolo fra p1 ed p2: %f\n", rad2deg(angolo2));

%% 4

% creiamo la funzione che applica il metodo di del gradiente coniugato e
% restituisce la matrice con tutte le soluzioni per ogni iterazione
% vedi file conjgrad_it.m

%% 5

% inizializzo guess iniziale, il numero di iterazioni massimo e la
% tolleranza sul residuo normalizzato come da testo (cambio la tolleranza
% come dalle soluzioni)
x0 = zeros(n, 1);
nmax = 1e3;
toll = 1e-14;

% chiamo la funzione creata al punto precedente
[x_gc, k_gc] = conjgrad_it(A, b, x0, nmax, toll);
fprintf("Iterazioni GC: %d\n", k_gc)

% verifico la correttezza del calcolo
% x = A \ b;

%% 6

% definissco la soluzione esatta
x_ex = ones(n, 1);

% inizializzo i vettori
err_rel_it = zeros(1, k_gc + 1);
r = zeros(n, k_gc + 1);
res_norm_it = zeros(1, k_gc + 1);
err_stima_it = zeros(1, k_gc + 1);

% calcolo per ogni iterazione l'errore relativo e il residuo normalizzato
for k = 0:k_gc
    err_rel_it(k+1) = norm(x_ex - x_gc(:, k+1)) / norm(x_ex);
    r(:, k+1) = b - A*x_gc(:, k+1);
    res_norm_it(k+1) = norm(r(:, k+1)) / norm(b);
    err_stima_it(k+1) = K_A * res_norm_it(k+1);
end

semilogy(0:k_gc, err_rel_it, '-o', 0:k_gc, res_norm_it, '-o', ...
    0:k_gc, err_stima_it, '-o')
grid on
legend("Errore relativo", "Residuo normalizzato", "Stima")

% lo stimatore dell'errore basato sul residuo normalizzato relativo
% sottostima l’errore relativo vero per la maggior parte delle iterate
% la qualità del criterio d'arresto è poco soddisfacente in quanto vale
% la stima dell’errore
% err_rel_stima <= K(A) * res_norm
% dato che K_2(A) = K(A) = 8.9687e4 >> 1, l'errore vero può essere
% sigificativamente sottostimato dal residuo normalizzato
% la stima è però sempre valida