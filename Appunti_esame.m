% Esame

% Matrice tridiagonale
A = diag(ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);

% Matrice pentadiagonale
A = diag(ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1)...
    + diag(ones(n-2,1),2) + diag(ones(n-2,1),-2);


%% Errori

% Errore assoluto
err = x - x_it

% Errore relativo
err_rel = norm(x - x_it) / norm(x)

% Errore in norma di A
err_norm_A = sqrt((x_it - x)' * A * (x_it - x))

% Stima errore relativo
stima_err_rel = res_norm * cond(A)

% Residuo
res = b - A * x_it

% Residuo relativo
res_rel = norm(res) / norm(b)

% Residuo normalizzato
res_norm = norm(b - A * x_it) / norm(b)


%% Fattori di abbattimento

% Numero di condizionamento
% Più il numero di condizionamento è alto più gli errori saranno alti,
% cioè il problema è mal posto
% Senza comando cond
condA = max(eig(A)) / min (eig(A))
condA = max(eig(P\A)) / min (eig(P\A))
condA = lambda/mu           % A SDP
condA = sqrt(lambda/mu)     % A non SDP

% Per d minori la convergenza è più veloce
d = (cond(A) - 1) / (cond(A) + 1)
d = (cond(P\A) - 1) / (cond(P\A) + 1)


% Fattore di abbattimento
% Calcolo del fattore di abbattimento dell'errore dopo k iterazioni 
% ||x - xk||_A <= d^k ||x - x0||_A;
% k parametro dato: numero di iterazioni
d = (cond(A)-1) / (cond(A)+1);
f_abb_k = d^k
%
% Numero minimo di iterazioni
k = log2(f_abb_k) / log2(d)

% Fattore di abbattimento precondizionati
% Calcolo del fattore di abbattimento dell'errore dopo k iterazioni 
% ||x - xk||_A <= d^k ||x - x0||_A;
% k parametro dato: numero di iterazioni
K = max(eig(P\A)) / min(eig(P\A));
d = (K-1) / (K+1);
f_abb_k = d^k

% Gradiente coniugato non precondizionato
% Calcolo del fattore di abbattimento dell'errore dopo k iterazioni 
% ||x-xk||_A <= (2*c^k) / (1 + c^(2k)) * ||x-x0||_A;
% k parametro dato: numero di iterazioni
c = (sqrt(cond(A)) - 1) / (sqrt(cond(A)) + 1);
f_abb_k = (2 * c^k) / (1 + c^(2*k))

% Gradiente coniugato precondizionato
% Calcolo del fattore di abbattimento dell'errore dopo k iterazioni 
% ||x-xk||_A <= (2*c^k) / (1 + c^(2k)) * ||x-x0||_A;
% k parametro dato: numero di iterazioni
K = max(eig(P\A)) / min(eig(P\A));
c = (sqrt(K) - 1) / (sqrt(K) + 1);
f_abb_k = (2 * c^k) / (1 + c^(2*k))

% Stima errore per metodo del gradiente:
err0 = sqrt ((x - x0)' * A * (x - x0));
stima_err_g = fatt_abb_k * err0;


%% Convergenza

% Convergenza Jacobi per ogni x0 iniziale:
% Precondizionatore per Jacobi è la matrice diagonale estratta da A
% Calcolo la matrice di iterazione B_j
% Calcolo il raggio spettrale (rho_j) della matrice di iterazione che è il
% massimo del valore assoluto degli autovalori della matrice di iterazione
% Se:
% rho_j < 1 -> converge
% rho_j > 1 -> non converge
D = diag(diag(A));
B_j = eye(n) - D\A;
rho_j = max(abs(eig(B_j)))

% Convergenza Gauss-Seidel per ogni x0 iniziale
% Precondizionatore per Gauss-Seidel è la matrice triangolare
% inferiore estratta da A
% Calcolo la matrice di iterazione B_gs
% Calcolo il raggio spettrale (rho_gs) della matrice di iterazione che è il
% massimo del valore assoluto degli autovalori della matrice di iterazione
% Se:
% rho_j < 1 -> converge
% rho_j > 1 -> non converge
T = tril(A);
B_gs = eye(n) - T\A;
rho_gs = max(abs(eig(B_gs)))

% Condizioni sufficienti per la convergenza di Jacobi e Gauss-Seidel
% - A è non singolare e a dominanza diagonale stretta per righe:
%   Gauss-Seidel e Jacobi convergono per ogni x0
% - A è SDP (simmetrica definita positiva):
%   Gauss-Seidel converge per ogni x0
% - A è non singolare e tridiagonale, allora rho(B_gs)=rho(B_j)^2: 
%   sono entrambi convergenti o divergenti
%   Se convergono (rho<1) allora Guass-Seidel converge più velocemente


% Newton
% E' necessario scegliere x(0) “sufficientemente” vicino allo zero alpha
% per essere certi che converga
% 
% m=1 -> p=2
% m>1 -> p=1, 0<mu<1

% Newton modificato
% Se x0 è sufficientemente vicino ad alpha allora
% m=2 -> p=2


% Metodo delle iterazioni di punto fisso:
% Se abs(f'(alpha)) < 1 e x0 è sufficientemente vicino ad alpha,
% il metodo converge ad alpha (almeno linearmente)
%
% Se f'(alpha)=0, f''(alpha)~=0 e x0 è sufficientemente vicino ad alpha,
% il metodo converge  quadraticamente ad alpha (con ordine 2)
%
% Convergenza monotona se f'(alpha)>0 (da unire alla soluzione della
% convergenza)
%
% p = m+1;


% Gradiente coniugato precondizionato
% Se A e P simmetriche definite positive converge per ogni x0
% d = (cond(P\A) - 1) / (cond(P\A) + 1)

% Convergenza gradiente coniugato non precondizionato:
% IN ARITMETICA ESATTA converge in al più n iterazioni per ogni x0


% Richardson stazionario:
%  
% Condizione necessaria e sufficiente:
% rho(B_alpha) < 1
% con B_alpha = eye(n) - alpha * P \ A
%     rho_B_alpha = max(abs(eig(B_alpha)));
% Quindi deve valere che
% abs(eig((eye(n) - alpha * P\A)) < 1
% 
% Matrice A non SDP: converge per ogni x0 se e solo se 
alpha < 2 * real(abs(eig(P\A))) / abs(eig(P\A))^2
%
% Matrici A e P simmetriche e definite positive
% il metodo converge per ogni x0 se e solo se
0 < alpha < 2 / max(eig(P\A))
%
% Matrici A e P non-singolari e con autovalori tutti a valori reali
% allora il metodo converge per ogni x0 se e solo se
0 < alpha * eig(P\A)< 2
%
% Alpha_opt:
alpha_opt = 2 / (min(eig(P\A)) + max(eig(P\A)))

% Richardson dinamico - gradiente
% il metodo converge se A è simmetrica definita positiva per ogni
% scelta di x0


% Metodo delle potenze
%per studiare l'ordine di convergenza p e il fattore di convergenza
%  asintotico c utilizzo la funzione
[p, c] = stimap(lambda_vec)
% Prende in input il vettore soluzione del metodo ad ogni iterazione
[p, c] = stimap(lambda_vec);
p = p(end);    % ordine di convergenza
c = c(end);    % fattore di convergenza
%
% Con autovalori complessi coniugati la convergenza col metodo delle
% potenza non è assicurata
% - mu=i -> converge all'autovalore con parte immaginaria positiva
% - mu=-i -> converge all'autovalore con parte immaginaria negativa


% Funzione
% Convergenza: abs(phi'(x)) < 1
% Unicità di alpha: phi(x) appartiene a [a,b]


% Convergenza per matrice a caso: al più n^2 operazioni


%% Bisezione

% Quante iterazioni > 0 sono necessarie al metodo di bisezione al
% fine di garantire un errore inferiore a tot?
nmax = ceil(log2((b - a) / toll) - 1)


% Applicabilità per alpha:
% Le radici di molteplicità pari sono tali che f non cambia segno in un
% intorno di alpha:
% per radici di molteplicità pari il metodo di bisezione non è applicabile


% Verificare che il metodo sia applicabile
% Devo controllare che ci sia almeno uno zero quindi dovrò avere un y
% positivo e un y negativo e il loro prodotto dovrà essere minore di zero
if (f(a) * f(b)) < 0
    disp('Si puo'' applicare il metodo di bisezione');
else 
    disp('Non si puo'' applicare il metodo di bisezione \n');
end


%% Cerchi di Gershgorin

% Gli autovalori si possono localizzare sul piano complesso di A
% attraverso i cerchi di Gershgorin
% - gli autovalori di A appartengono alla regione di intersezione dei
%   cerchi colonna e i cerchi riga
% - se due regioni di spazio S1 e S2 sono disgiunte, allora S1 contiene
%   esattamente m autovalori e S2 n-m autovalori
% - se un cerchio è isolato, contiene un solo autovalore
% - raggio = somma elementi sulla riga/colonna escluso l'elemento sulla
%   diagonale
% - elemento sulla diagonale = centro del cerchio


%% FATTORIZZAZIONE LU:

% Scompone A come A = L*U, con L matrice triangolare inferiore e U
% matrice trangolare superiore
% A * x = b
% L * U * x = b
% Ottengo i sistemi:
% L * y = b
% U * x = y


% Condizione necessaria sufficiente:
% A invertibile è LU fattorizzabile sse tutti i minori principali di A sono
% invertibili
% Tutte le sottomatrici principali di A devono essere non singolari,
% cioè il loro determinante deve essere diverso da 0
dA = 1;
% nA = numero di colonne
nA = size(A, 1);
% Ciclo while da cui esco se ho analizzato tutti gli elementi o se trovo
% che il determinante di una sottomatrice è uguale a 0
while ( (i < nA) && (dA ~= 0) )
    dA = det( A(1:i, 1:i) );
    i = i + 1;        
end
% Se dA è uguale a 0, vuol dire che sono uscito dal ciclo perché ho trovato
% un determinante uguale a 0 e di conseguenza la matrice non è LU
% fattorizzabile
if (dA == 0)
    disp('Non posso applicare la fattorizzazione LU')
% Se dA è diverso da 0, vuol dire che ho analizzato tutti i determinanti e
% nessuno è uguale a 0. Di conseguenza posso applicare la fattorizzazione
% LU alla matrice
else
    disp('Posso applicare la fattorizzazione LU')
end


% Condizioni sufficienti per la fattorizzabilità data A invertibile:
%    - A è a dominanza stretta per righe
%    - A è dominanza stretta per colonne
%    - A è SDP (simmetrica definita positiva)

% Verifico se A è a dominanza diagonale stretta per righe:
% Trovo la diagonale di A e la salvo in un vettore colonna
A_diag = diag(abs(A));
% Trovo la somma degli elementi per riga esclusa la diagonale e li salvo
% in un vettore colonna
A_outdiag = sum(abs(A),2) - A_diag;
% Se le somme sono minori del corrispondente elemento sulla diagonale,
% la matrice è a dominanza diagonale per righe
if (A_diag > A_outdiag)
    disp('A è a dominanza diagonale per righe')
end

% Verifico se A è a dominanza diagonale stretta per colonne:
% Trovo la diagonale di A e la salvo in un vettore riga
A_diag = diag(abs(A))';
% Trovo la somma degli elementi per colonna esclusa la diagonale e li salvo
% in un vettore riga
A_outdiag = sum(abs(A),1) - A_diag;
% Se le somme sono minori del corrispondente elemento sulla diagonale,
% la matrice è a dominanza diagonale per colonne
if (A_diag > A_outdiag)
    disp('A è a dominanza diagonale per colonne')
end

% Verifico se A è SDP (simmetrica definita positiva):
% Verifico che A sia simmetrica, ossia che A' (A trasposta) sia
% uguale ad A:
if(A == A')
    % Se A è simmetrica, allora verifico che sia definita positiva,
    % ossia che i suoi autovalori siano positivi:
    if(eig(A) > 0)
        % Se entrambe le condizioni sono rispettate visualizzo il
        % seguente messaggio:
        disp('Matrice è SDP')
    end
end


%% Soluzione sistemi

% Fattorizzazione LU
y = L \ b;
x = U \ y;

[L, U, P] = lu(A)
y = fwsub (L, P * b);
x = bksub (U, y);

% Pivoting
P * A * x = P * b
U * x = y

L * y = P * b
U * x = y

% MEG con pivoting totale
[L, U, P, Q] = lu(sparse(A));
L * y_star = P * b
U * x_star = y_star
x = Q * x_star


% Cholesky
% Trovo R matrice della fattorizzazione di Cholesky tale che
% A = R' * R
R = chol(A);
% Risoluzione del sistema:
% R' * y = b
% R * x = y
y = R' \ b;
x = R \ y;


% Fattorizzazione QR
% A -> m*m, x -> n*1, b -> m*1
% Uso il comando qr per trovare Q e R tali che A = Q * R
[Q, R] = qr(A);
% Q -> m*m
% R -> m*n matrice rettangolare i cui elementi sotto la diagonale
% principale sono nulli
% Risolvo il sistema
% R * x = Q' * b;
x = R \ (Q' * b);


% Gradiente coniugato non precondizionato
[x_pcg,~,~,iter] = pcg(A,b,tol,nmax,eye(n),eye(n),x0);
[xk, it] = conjgrad_it(A, b, x0, nmax, toll);
% Precondizionato
[x_pcg,~,~,iter] = pcg(A,b,tol,nmax,P,eye(n),x0);


% Metodo iterativo lineare
B = eye(n) - P \ A;
g = P \ b;


%% Grafici

plot(x,'bo-', 'LineWidth', 2)

semilogy(N, It_np, N, It_p, 'Linewidth', 2);

% Per fare griglia
grid on;

% Per fare legenda
legend(' ', ' ');

% Asse x
y0 = zeros(spaziatura, 1);         
hold on
plot(x, y0)

% Per fargli fare tre grafici diversi
figure(1)
spy(A)
figure(2)
spy(L)
figure(3)
spy(U)


%%  Teoria

% F = F(b, t, L, U)
%   - b = base
%   - t = numero di cifre significative
%   - L = esponente minimo
%   - U = esponente massimo

% Epsilon macchina = b^(1-t)

% Errore di arrotondamento:
% ||x - fl(x)|| / ||x|| = 1/2 * epsilon macchina

% x = (-1)^s * m * b^(e-t)

%   - s = segno
%   - m = mantissa
%   - e = esponente (L <= e <= U)

% Numero massimo rappresentabile
% n = m * b^(U-t)

% Numero massimo
% x_max = b^U * (1 - b^(-t))

% Numero minimo
% x_min = b^(L-1)

% Cardinalità
% card(F0) = 2 * (b-1) * b^(t-1) * (U-L+1)


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

%% Altro

% Se A è invertibile -> L e U sono invertibili
% det(A) = det(L*U) = det(L)*det(U) = det(U)

% Conversione da binario a decimale:
bin2dec(' ')


%% Da ricordare

% Cambiare nome alla x data in output da una funzione