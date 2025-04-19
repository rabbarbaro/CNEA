%% matrici tri/penta diagonali

% Matrice tridiagonale
A = diag(ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);

% Matrice pentadiagonale
A = diag(ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1)...
    + diag(ones(n-2,1),2) + diag(ones(n-2,1),-2);

%% aritmetica floating point

% F = F(b, t, L, U)
%   - b = base
%   - t = numero di cifre significative
%   - L = esponente minimo
%   - U = esponente massimo

%   - s = segno
%   - m = mantissa
%   - e = esponente (L <= e <= U)

% numero in floating point, ricorda mantissa in decimale con bin2dec()
x = (-1)^s * m * b^(e-t);
% epsilon macchina
eps_M = b^(1 - t);
% errore arrotondamento
abs(x - fl(x)) / abs(x) <= 1/2 * eps_M;
% numero minimo floating point
x_min = b^(L - 1);
% numero massimo floating point
x_max = b^U * (1 - b^(-t));
% cardinalità insieme (numeri positivi e negativi, escluso lo 0)
card = 2 * (b-1) * b^(t-1) * (U-L+1);

%% metodi diretti

% sistema diagonale
%   x_ii = b_ii / d_ii

% sistemi triangolari
y = fwsub(L, b);
x = bksub(U, y);

% fattorizzazione LU
% Scompone A come A = L*U, con L matrice triangolare inferiore e U
% matrice trangolare superiore
%   A * x = b
%   L * U * x = b
% Ottengo i sistemi:
%   L * y = b
%   U * x = y
% applicando il pivoting
%   P * A * x = P * b
%   L * U * x = P *b
%   L * y = P * b
%   U * x = y
[L, U] = lu(A);
y = L \ b;
x = U \ y;
% con pivoting
[L, U, P] = lu(A);
y = fwsub (L, P * b);
x = bksub (U, y);
% con pivoting totale
[L, U, P, Q] = lu(sparse(A));
L * y_star = P * b;
U * x_star = y_star;
x = Q * x_star;

% % condizione necessaria e sufficiente per fattorizzazione LU:
% % A invertibile è LU fattorizzabile sse tutti i minori principali di A sono
% % invertibili (tutte le sottomatrici principali di A devono essere non
% % singolari, cioè il loro determinante deve essere diverso da 0)
% dA = 1;
% % nA = numero di colonne
% nA = size(A, 1);
% % Ciclo while da cui esco se ho analizzato tutti gli elementi o se trovo
% % che il determinante di una sottomatrice è uguale a 0
% while ( (i < nA) && (dA ~= 0) )
%     dA = det( A(1:i, 1:i) );
%     i = i + 1;        
% end
% % Se dA è uguale a 0, vuol dire che sono uscito dal ciclo perché ho trovato
% % un determinante uguale a 0 e di conseguenza la matrice non è LU
% % fattorizzabile
% if (dA == 0)
%     disp('Non posso applicare la fattorizzazione LU')
% % Se dA è diverso da 0, vuol dire che ho analizzato tutti i determinanti e
% % nessuno è uguale a 0. Di conseguenza posso applicare la fattorizzazione
% % LU alla matrice
% else
%     disp('Posso applicare la fattorizzazione LU')
% end

% % Condizioni sufficienti per fattorizzazione LU (A invertibile):
% %    - A è a dominanza stretta per righe
% %    - A è dominanza stretta per colonne
% %    - A è SDP (simmetrica definita positiva)
% 
% % Verifico se A è a dominanza diagonale stretta per righe:
% % Trovo la diagonale di A e la salvo in un vettore colonna
% A_diag = diag(abs(A));
% % Trovo la somma degli elementi per riga esclusa la diagonale e li salvo
% % in un vettore colonna
% A_outdiag = sum(abs(A),2) - A_diag;
% % Se le somme sono minori del corrispondente elemento sulla diagonale,
% % la matrice è a dominanza diagonale per righe
% if (A_diag > A_outdiag)
%     disp('A è a dominanza diagonale per righe')
% end
% 
% % Verifico se A è a dominanza diagonale stretta per colonne:
% % Trovo la diagonale di A e la salvo in un vettore riga
% A_diag = diag(abs(A))';
% % Trovo la somma degli elementi per colonna esclusa la diagonale e li salvo
% % in un vettore riga
% A_outdiag = sum(abs(A),1) - A_diag;
% % Se le somme sono minori del corrispondente elemento sulla diagonale,
% % la matrice è a dominanza diagonale per colonne
% if (A_diag > A_outdiag)
%     disp('A è a dominanza diagonale per colonne')
% end
% 
% % Verifico se A è SDP (simmetrica definita positiva):
% % Verifico che A sia simmetrica, ossia che A' (A trasposta) sia
% % uguale ad A:
% if(A == A')
%     % Se A è simmetrica, allora verifico che sia definita positiva,
%     % ossia che i suoi autovalori siano positivi:
%     if(eig(A) > 0)
%         % Se entrambe le condizioni sono rispettate visualizzo il
%         % seguente messaggio:
%         disp('Matrice è SDP')
%     end
% end

% thomas
% se la matrice è tridiagonale
[L, U, x] = thomas(A, b);

% cholesky
% se la matrice è SDP
% Trovo R matrice della fattorizzazione di Cholesky tale che A = R' * R
R = chol(A);
%   R' * y = b
%   R * x = y
y = R' \ b;
x = R \ y;

%% sistemi sovradeterminati

% fattorizzazione QR
%   - A -> m*n
%   - x -> n*1
%   - b -> m*1
% trovo soluzione ai minimi quadrati
%   A' * A * x_star = A' * b
% Uso il comando qr per trovare Q e R tali che A = Q * R
[Q, R] = qr(A);
% Q -> m*m
% R -> m*n matrice rettangolare i cui elementi sotto la diagonale
%           principale sono nulli
% Risolvo il sistema
%   R * x_star = Q' * b;
x_star = R \ (Q' * b);      % bksub

%% metodi iterativi

% Metodo iterativo
%   x_k+1 = B * x_k + g
B = eye(n) - P \ A;
g = P \ b;

% convergenza metodo iterativo
%   ||e_k|| = ||B^k * e_0|| <= ||B^k||*||e_0|| <= ||B||^k*||e_0||
%   ||e_k|| <= (rho(B))^k*||e_0||
% converge se rho(B) < 1

% residuo precondizionato
%   P * z = r
% metodo iterativo
% dato x0, r0 = b - A * x0
% while
%   P * z_k = r_k
%   x_k+1 = x_k + z_k
%   r_k+1 = r_k - A * z_k

% jacobi
% Precondizionatore per Jacobi è la matrice diagonale estratta da A
% Calcolo la matrice di iterazione B_j
% Calcolo il raggio spettrale (rho_j) della matrice di iterazione che è il
% massimo del valore assoluto degli autovalori della matrice di iterazione
% Se:
%   rho_j < 1 -> converge
%   rho_j > 1 -> non converge
D = diag(diag(A));
B_j = eye(n) - D\A;
rho_j = max(abs(eig(B_j)));

% gauss-seidel
% Precondizionatore per Gauss-Seidel è la matrice triangolare
% inferiore estratta da A
% Calcolo la matrice di iterazione B_gs
% Calcolo il raggio spettrale (rho_gs) della matrice di iterazione che è il
% massimo del valore assoluto degli autovalori della matrice di iterazione
% Se:
%   rho_j < 1 -> converge
%   rho_j > 1 -> non converge
T = tril(A);
B_gs = eye(n) - T\A;
rho_gs = max(abs(eig(B_gs)));

% Condizioni sufficienti per la convergenza di Jacobi e Gauss-Seidel
%   - A è non singolare e a dominanza diagonale stretta per righe:
%       Gauss-Seidel e Jacobi convergono per ogni x0
%   - A è SDP (simmetrica definita positiva):
%       Gauss-Seidel converge per ogni x0
%   - A è non singolare e tridiagonale, allora rho(B_gs)=rho(B_j)^2: 
%       sono entrambi convergenti o divergenti, se convergono (rho<1)
%       allora Guass-Seidel converge più velocemente

% richardson stazionario:
B_alpha = eye(n) - alpha * P \ A;
g_alpha = alpha * P \ b;
% Condizione necessaria e sufficiente:
%   rho(B_alpha) < 1
%   rho_B_alpha = max(abs(eig(B_alpha)));
% Quindi deve valere che
%   abs(eig((eye(n) - alpha * P\A)) < 1
% Matrice A non SDP: converge per ogni x0 se e solo se:
alpha < 2 * real(abs(eig(P\A))) / abs(eig(P\A))^2;      % per ogni autovalore
% Matrici A e P non-singolari e con autovalori tutti a valori reali
% allora il metodo converge per ogni x0 se e solo se:
0 < alpha * eig(P\A) < 2;       % per ogni autovalore
% Matrici A e P SDP: converge per ogni x0 se e solo se:
0 < alpha < 2 / max(eig(P\A))
% Alpha_opt:
alpha_opt = 2 / (min(eig(P\A)) + max(eig(P\A)));

% richardson dinamico - gradiente:
alpha_k = (z' * r) / (z' * A * z);
% se non è precondizionato z = r
% il metodo converge se A è simmetrica definita positiva per ogni
% scelta di x0

% gradiente coniugato
% definisco anche il parametro beta_k durante il metodo (non richardson),
% ottimalità garantita rispetto a tutte le direzioni precedenti
% Non precondizionato
[x_pcg,~,~,iter] = pcg(A,b,tol,nmax,eye(n),eye(n),x0);
[xk, it] = conjgrad_it(A, b, x0, nmax, toll);
% Precondizionato
[x_pcg,~,~,iter] = pcg(A,b,tol,nmax,P,eye(n),x0);
% Convergenza gradiente coniugato non precondizionato:
% IN ARITMETICA ESATTA converge in al più n iterazioni per ogni x0

%% convergenza richardson (gradiente) e gradiente coniugato (fattore di abbattimento)

% richardson (gradiente):
%   ||x - xk||_A <= d^k ||x - x0||_A;
% k parametro dato: numero di iterazioni
d = (cond(A)-1) / (cond(A)+1);
f_abb_k = d^k;
% Numero minimo di iterazioni data tolleranza
k = log(tol) / log(d);
% Fattore di abbattimento precondizionati
K_invPA = max(eig(P\A)) / min(eig(P\A));
d = (K-1) / (K+1);
f_abb_k = d^k;

% Gradiente coniugato:
% non precondizionato
%   ||x-xk||_A <= (2*c^k) / (1 + c^(2k)) * ||x-x0||_A;
% k parametro dato: numero di iterazioni
c = (sqrt(cond(A)) - 1) / (sqrt(cond(A)) + 1);
f_abb_k = (2 * c^k) / (1 + c^(2*k));
% Gradiente coniugato precondizionato
K_invPA = max(eig(P\A)) / min(eig(P\A));
c = (sqrt(K) - 1) / (sqrt(K) + 1);
f_abb_k = (2 * c^k) / (1 + c^(2*k));
% per k grandi vale circa
f_abb_k = 2 * c^k;

% Stima errore:
err0_normA = sqrt ((x - x0)' * A * (x - x0));
stima_err_normA = fatt_abb_k * err0_normA;

%% numeri di condizionamento, raggio spettrale

% raggio spettrale
rho_A = max(abs(eig(A)));
% numero di condizionamento (2)
K2_A = cond(A);
K2_A = sqrt(max(abs/eig(A'*A))/min(abs(eig(A'*A))));
% numero di condizionamento spettrale
K_A = rho_A * rho_invA;
K_A = max(abs(eig(A)))/min(abs(eig(A)));
% se A SDP
K2_A = K_A = max(eig(A))/min((eig(A));

%% errore e residuo assoluto e relativo

% Errore assoluto
err = x_ex - x_it;
% Errore relativo
err_rel = norm(x_ex - x_it) / norm(x_ex);
% Errore in norma di A
err_norm_A = sqrt((x_it - x_ex)' * A * (x_it - x_ex))
% Residuo
res = b - A * x_it;
% Residuo relativo/normalizzato
res_rel = norm(res) / norm(b);
res_norm = norm(b - A * x_it) / norm(b);

% Stima (maggiorazione) errore relativo (quando deltaA ~ 0, wilkinson-higham)
stima_err_rel = res_norm * cond(A);

% per i metodi iterativi vale
stima_err_rel_k = res_norm_K * cond(A);
norm(err_k) = 1/(1-rho_B) * norm(delta_k);

%% autovalori e autovettori

% quoziente di rayleigh
q = (x' * A * x) / (x' * x);
% metodo delle potenze
% trova autovalore di massimo modulo se ben distinto
[lambda,x,iter] = eigpower(A,tol,nmax,x0);
% metodo delle potenze inverse
% trova autovalore di minimo modulo se ben distinto
[lambda,x,iter] = invpower(A,tol,nmax,x0);
% ATTENZIONE SE SI SCEGLIE UN AUTOVETTORE DI PARTENZA PERPENDICOLARE A
% QUELLO CORRISPONDENTE ALL'AUTOVETTORE DA TROVARE
% metodo delle potenze inverse con shift
% trova autovalore più vicino allo shift se ben distinto
[lambda,x,iter] = invpowershift(A,mu,tol,nmax,x0);

% i metodi delle potenze convergono sempre con ordine 1
% matrice A generica
%   e_k = |lambdaA(1) - lambda_k| <= C * |(eigA(2) - /eigA(1)|^k
% matrice A hermitiana (simmetrica)
%   e_k = |lambdaA(1) - lambda_k| <= C * |(eigA(2) - /eigA(1)|^2k

% se viene chiesto di trovare il numero minimo di iterazioni
%   |lambdaA(1) - lambda_k| / |lambdaA(1) - lambda_0| < toll
% dato che se k = 0 vale:
%   |lambdaA(1) - lambda_0| <= C
% diventa (se hermitiana)
%   |(eigA(2)/eigA(1)|^2k < tol
%   2k * log(|(eigA(2)/eigA(1)|) < log(tol)
%   dato che log(|(eigA(2)/eigA(1)|) < 0 cambio verso alla disequazione
%   N = 0.5 * log(tol) / log(abs(eigA(2)/eigA(1)));

% Cerchi di Gershgorin
% Gli autovalori si possono localizzare sul piano complesso di A
% attraverso i cerchi di Gershgorin
%   - gli autovalori di A appartengono alla regione di intersezione dei
%       cerchi colonna e i cerchi riga
%   - se due regioni di spazio S1 e S2 sono disgiunte, allora S1 contiene
%       esattamente m autovalori e S2 n-m autovalori
%   - se un cerchio è isolato, contiene un solo autovalore
%   - raggio = somma elementi sulla riga/colonna escluso l'elemento sulla
%       diagonale
%   - elemento sulla diagonale = centro del cerchio

% iterazioni qr
D = qrbasic(A,tol,nmax);
%   err_est_k = max_(i=2:n, j=1:i-1)(A_tilde(i, j)
% ovvero il massimo dei pivot

%% SVD

[U, S, V] = svd(A);
% per il rapporto di compressione, scelto k per l'approssimazione a rango ridotto
[m,n] = size(A);
uncompressed_size = m*n;
A_k = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
compressed_size = k * (m+n+1);

% vedi teorema di schidt-eckart-young per norma2 e normaFro (cattura dell'
% "energia" della matrice

%% equazioni non lineari

% ordine di convergenza p
%   lim_k->inf |x_k+1 - alpha|/|x_k - alpha|^p = mu
% se p = 1, mu in (0,1)

% bisezione
[xvect, it] = bisez(a,b,toll,fun);
% numero di iterazioni necessarie per garantire un errore inferiore a tol
nmax = ceil(log2((b - a) / toll) - 1);
% Applicabilità per alpha:
%  - Le radici di molteplicità pari sono tali che f non cambia segno in un
%       intorno di alpha:
%  - per radici di molteplicità pari il metodo di bisezione non è applicabile
% Verificare che il metodo sia applicabile
% Devo controllare che ci sia almeno uno zero quindi dovrò avere un y
% positivo e un y negativo e il loro prodotto dovrà essere minore di zero
if (f(a) * f(b)) < 0
    disp('Si puo'' applicare il metodo di bisezione');
else 
    disp('Non si puo'' applicare il metodo di bisezione \n');
end

% Newton
x_k+1 = x_k - f(x_k)/df(x_k);
% Bisogna scegliere x(0) "sufficientemente" vicino allo zero alpha per
% essere certi che converga
% per zero semplice (m=1 -> p=2)
%   lim_k->inf (x_k+1 - alpha)/(x_k - alpha)^2 = 1/2 * d2f(alpha)/df(alpha) = mu
% per zero multiplo (m>1 -> p=1)
%   lim_k->inf (x_k+1 - alpha)/(x_k - alpha) = mu

% newton modificato
x_k+1 = x_k - m * f(x_k)/df(x_k);
% sempre p = 2

% newton adattivo
m = (x_k-1 - x_k-2) / (2*x_k-1 - x_k - x_k-2);

% quasi-newton
x_k+1 = x_k - f(x_k)/q_k;
% newton
q_k = df(x_k);
% newton modificato
q_k = df(x_K)/m;
% metodo delle corde (p = 1 se zero semplice con condizioni su q_c, non garantito per zeri multipli)
q_k = q_C = (f(b) - f(a)) / (b - a);
% metodo delle secanti (p = (1 + sqrt(5))/2 se zero semplice, p = 1 se zero multiplo
q_k = (f(x_k) - f(x_k-1)) / (x_k - x_k-1);
%   q0 per secanti scelto come corde

%% newton per sistemi di equazioni non lineari

% % metodo completo
% % definisco la funzione vettoriale e la sua jacobiana
% F = @(x1, x2) ...
% JF = @(x1, x2) ...
% % iterata iniziale
% x = [...]';
% % numero massimo di iterazioni
% kmax = ...;
% % inizializzo xvec come vettore e l'iterata 1
% xvec = [];
% % per ogni iterazione
% % bisognerebbe anche verificare a ogni iterazione che Jmat(x) corrente
% sia diversa da 0 per ogni x
% for k = 1:kmax
%     % valuto funzione e jacobiana nell'iterata corrente
%     Fvec = F(x(1), x(2));
%     Jmat = JF(x(1), x(2));
%     % risolvo il sistema lineare
%     x = x - Jmat\Fvec;
%     % salvo l'iterata
%     xvec = [xvec x];
% end
% % soluzione trovata
% xvec(:, end);
% % verifico se il determinante della jacobiana sia != 0
% det(JF(alpha(1), alpha(2)));
% % se è != 0 ci aspettiamo convergenza quadratica

% % per stimare l'ordine di convergenza di newton per sitemi
% % calcoliamo il vettore degli errori
% errvec = [];
% for k = 1:kmax-1
%     errvec = [errvec, norm(xvec(:, k) - alpha)];
% end
% % calcoliamo un vettore contente il rapporto fra errori successivi
% % elevando il denominatore all'ordine di convergenza troveremo il fattore
% % di convergenza
% % aumentando l'ordine di convergenza, il primo che avrà un fattore di
% % convergenza non nullo o non infinito sarà l'ordine di convergenza
% % vale:
% % lim ||x^(k+1)-alpha||/||x^(k)-alpha||^p = mu
% % p = 1?
% conv_ord1 = errvec(2:end) ./ errvec(1:end-1)
% % se la convergenza di ordine 1 tende a 0 --> ordine superiore
% % p = 2?
% conv_ord2 = errvec(2:end) ./ errvec(1:end-1).^2
% % e così via

%% iterazioni di punto fisso

% algoritmo di punto fisso
x_k+1 = phi(x_k);
% phi più semplice possibile
phi = @(x) f(x) + x;
% se converge il limite della succesione è un punto fisso

% convergenza globale in un intervallo [a, b]
%   1. se phi in C0, immagine in [a, b], allora esiste almeno un punto
%       fisso alpha in [a, b]
%   2. se esiste costante L in [0, 1) t.c. |phi(x1) - phi(x1)| <= L * |x1 - x2|
%       per ogni x1, x2 in [a, b] allora alpha è unico in [a, b] e il
%       metodo delle iterazioni di punto fisso converge per ogni x0
err_k <= L^k * err_0;
% convergenza globale in un intervallo [a, b]
% se phi in C1, immagine in [a, b], |dphi(x)| < 1 per ogni x1, x2 in [a, b]
% allora alpha è unico in [a, b] e il metodo delle iterazioni di punto
% fisso converge per ogni x0 con ordine almeno p = 1
%   lim_k->inf (x_k+1 - alpha)/(x_k - alpha) = dphi(alpha)
% dove dphi(alpha) è il fattore di convergenza asintotico

% convergenza locale intorno punto fisso (1)
% se phi in C1, I_alpha intorno di alpha, |dphi(alpha)| < 1 se x0 è
% "sufficientemente" vicino ad alpha allora il metodo delle iterazioni di
% punto fisso converge con ordine almeno p = 1
%   lim_k->inf (x_k+1 - alpha)/(x_k - alpha) = dphi(alpha)
% dove dphi(alpha) è il fattore di convergenza asintotico
%   - se |dphi(alpha)| < 1 converge
%   - se |dphi(alpha)| = 1 dipende
%   - se |dphi(alpha)| > 1 non converge (a meno che x0 = alpha)

% convergenza locale intorno punto fisso (2)
% se phi in Cp, I_alpha intorno di alpha, |d(i)phi(alpha)| = 0 per ogni i =
% 1:p-1, |d(p)phi(alpha)| != 0, se x0 è "sufficientemente" vicino ad alpha
% allora il metodo delle iterazioni di punto fisso converge con ordine p
%   lim_k->inf (x_k+1 - alpha)/(x_k - alpha)^p = 1/p! * d(p)phi(alpha)
% dove 1/p! * d(p)phi(alpha) è il fattore di convergenza asintotico

% convergenza monotona se dphi(alpha)>0 (da mettere a sistema con la
% soluzione della convergenza)

% usando lagrange
%   x_k - alpha = - 1 / (1 - dphi(xi_k)) * delta_k
% lo stimatore d'errore basato su iterate successive è soddisfacente se
%   dphi(alpha) ~ 0
%   err_k ~ err_tilde_k+1
% lo stimatore d'errore basato su iterate successive è soddisfacente se
%   dphi(alpha) ~ -1
%   err_k ~ 1/2 * err_tilde_k+1
% lo stimatore d'errore basato su iterate successive è INSODDISFACENTE se
%   dphi(alpha) ~ 1
%   err_k << err_tilde_k+1

% newton come punto fisso
phi_N = @(x) x - f(x)/df(x);
% se f in Cm, I_alpha intorno di alpha, m >= 1 molteplicità di alpha, phi_N
% è tale che
dphi_N(alpha) = 1 - 1/m;

% se f in C2, alpha zero semplice (m = 1), x0 "sufficientemente" vicino ad
% alpha, newton converge quadraticamente (p = 2)
%   lim_k->inf (x_k+1 - alpha)/(x_k - alpha)^2 = 1/2 * d2phi_N(alpha) = 1/2 * d2f(alpha)/df(alpha)
%   mu = 1/2 * d2f(alpha)/df(alpha)

% se f in Cm, alpha zero multiplo (m > 1), x0 "sufficientemente" vicino ad
% alpha, newton converge linearmente (p = 1)
%   lim_k->inf (x_k+1 - alpha)/(x_k - alpha) = dphi_N(alpha) = 1 - 1/m != 0
%   mu = dphi_N(alpha) = 1 - 1/m

% se f in Cm, alpha zero (m >= 1), x0 "sufficientemente" vicino ad alpha,
% newton modificato converge quadraticamente (p = 2)
%   dphi_N(alpha) = 1 - m * 1/m = 0 sempre

% lo stimatore d'errore basato su iterate successive è soddisfacente se lo
% zero è semplice, INSODDISFACENTE se è multiplo (newton)
%   err_k = m*err_tilde_k+1 = m * |delta_k)
% lo stimatore d'errore basato su iterate successive è sempre soddisfacente
% (newton modificato)

% % se è richiesto m e non si vogliono fare le derivate
% % chiamo newton
% [xvect_N, it_N] = newton(x0, nmax, toll, f, df, 1);
% % stimo p e c
% [p_N, c_N] = stimap(xvect_N);
% % dalla teoria
% %   lim ||x^(k+1)-alpha||/||x^(k)-alpha||^p = 1/p! * d(p)phi(alpha)
% % se p = 1
% %   lim ||x^(k+1)-alpha||/||x^(k)-alpha|| = dphi(alpha)
% % dalla teoria
% %   dphi_newton(alpha) = 1 - 1/m
% %   m = 1 / (1 - dphi_newton(alpha))
% % dove dphi_newton(alpha) è proprio c_N, ovvero mu (fattore di convergenza asintotico)
% m = 1/(1 - c_N(end));
% % è possibile applicare newton modificato con m appena trovato
% [xvect_Nm, it_Nm] = newton(x0, nmax, toll, f, df, m);

% metodo delle corde come iterazioni di punto fisso
phi_C = @(x) x - f(x)/q_C;
q_C = (f(b) - f(a)) / (b - a);
%   dphi_C(x) = 1 - df(x)/q_C
% q_C > 1/2 * df(alpha)   -   se df(alpha) > 0
% q_C < 1/2 * df(alpha)   -   se df(alpha) < 0
% se zero multiplo, dphi_C(alpha) = 1, conevrgenza non garantita
% se |dphi_C(alpha)| < 1 convergenza lineare, mu = 1 - df(alpha)/q_C
% se dphi_C(alpha) = 0 (ovvero q_C = df(alpha)) convergenza quadratoca, mu = -d2f(alpha)/(2 * q_C)

% iterazioni di punto fisso per funzioni vettoriali
% come quello scalare ma con vettori, differenza nella convergenza
% lim_k->inf ||x_k+1 - alpha||/||x_k - alpha|| = rho(J_phi(alpha))

%% ottimizzazione numerica

% gradiente ed hessiana
% punto stazionario se gradiente nullo
% se hessiana definita positiva in punto stazionario è punto di minimo
% locale

% sezione aurea
phi = (1 + sqrt(5))/2;
[xv, k, errv] = sezione_aurea(f, a, b, tol, nmax);
% errore 
%   err_k <= err_tilde_k = |I_k|/2 = (b - a) / (2 * phi^k)

% discesa
x_k+1 = x_k + alpha_k * d_k;

% condizioni di wolfe
%   1. Phi(x_k + alpha_k * d_k) <= phi(x_k) + c1 * alpha_k * gradPhi(x_k)' * d_k
%   2. |gradPhi(x_k + alpha_k * d_k)' * d_k| <= c2 * |gradPhi(x_k)' * d_k|
% la 1. impone che la riduzioni di di Phi sia significativa rispetto alla
% variazione predetta dal gradiente (impedisce passi troppo grandi)
% la 2. assicura che il gradiente aggiornato k+1 non cambi troppo
% velocemente rispetto al gradiente k
% backtracking
[alpha_k, it] = backtracking(Phi, GradPhi, x_k, d_k, c1, rho, nmax);

% gradiente
d_k = -gradPhi(x_k);
% gradiente coniugato
d_0 = -gradPhi(x_0);
d_k+1 = -gradPhi(x_k+1) + beta_k * d_k;
% beta_k scelta con vari criteri
% gradiente e gradiente coniugato
[xvect, it] = gradiente_coniugato_opt(Phi, GradPhi, flag_metodo, x0, toll, nmax);
% gradiente converge con ordine p <= 1;
% gradiente coniugato converge con ordine p >= 1;

% newton
[xvect, it] = newton_opt(GradPhi, HessPhi, x0, toll, nmax);
% d_k = -(H_phi(x_k))^-1 * gradPhi(x_k)
% se det(H_phi(x_ex)) != 0 converge quadraticamente:
%   lim_k->inf ||x_k+1 - x_ex||/||x_k - x_ex||^2 = mu
% se det(H_phi(x_ex)) = 0 converge linearmente

% bfgs
% approssimazione hessiana con B
d_k = -B_k * gradPhi(x_k);
% con B_k definita nel metodo
[xvect, it] = bfgs(Phi, GradPhi, x0, toll, nmax);
% converge in modo superlineare

%% Scelta metodo

% Sistemi lineari
%   - MEG (con pivoting): per ogni A invertibile
%   - Thomas: A tridiagonale invertibile
%   - Cholesky: A SDP

% Fattorizzazione QR: sistemi sovradeterminati (rettangolari)

% Decomposizione additiva
%   - Jacobi
%   - Gauss-Seidel

% Richardson: A SDP (simmetrica definita positiva)
%   - Gradiente:
% metodi di richardson dinamico le direzioni di discesa scelte ad ogni
% iterazione per arrivare alla soluzione sono perpendicolari tra loro
%   - Gradiente coniugato: 
% le direzioni di discesa scelte ad ogni iterazione per arrivare alla
% soluzione non sono perpendicolari tra loro. La nuova iterata al passo
% k+1 deve essere ottimale non solo alla direzione precedente k,
% ma a tutte le direzioni usate fino a quel momento.
%
% - Stazionario: alpha = cost
% - Dinamico: alpha cambia ad ogni iterazione

% Autovalori e autovettori:
%   - Potenze dirette: cerca l'autovalore di modulo max
%             i due autovalori di ordine massimo devono essere distinti
%             in modulo
%   - Potenze inverse: cerca l'autovalore di modulo minimo
%             lambda_n distinto in modulo da lambda_n-1
%   - Potenze inverse con shift: cerca l'autovalore più vicino allo shift
%   - Metodo delle iterazioni qr: calcolare tutti gli autovalori
%             autovalori tutti reali e distinti in modulo

% Newton
% Il metodo di Newton è applicabile a una funzione f ∈ C0(I)
% differenziabile in I

%% Costi

% - Fattorizzazione LU: 2/3*n^3

% Per matrici triangolari
% A * x = y : costo = n^2

% - Thomas: 8*n-7
%   Di cui
%   3*(n-1) per scomposizione LU,
%   2*(n-1) per risolvere L*y=b
%   3*(n-1)+1 per risolvere U*x=y

% - Cholesky: 1/3*n^3 + 2*n^2
%   Di cui
%   1/3*n^3 per risolvere R = chol(A)
%   n^2 per risolvere R'y=b
%   n^2 per risolvere Rx=y

% - Fattorizzazione QR: m*n

% - Risolvere L (fwsub): n^2
% - Risolvere U (bksub): n^2

% Matrice sparsa: calcolare a mano