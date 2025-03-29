clear
clc
close all

A = [50 1 3
    1 6 0
    3 0 1];
B = [50 1 10
    3 20 1
    10 4 70];
C = [0 0 9
    5 4 3
    1 2 6];

%% 1. verifica condizioni per fattorizzazione LU

% questa cosa è abbastanza noiosa e pallosa, negli esercizi d'esame
% si potrebbe fare nel terminale

%% 1.1 A SDP?

% verifico se è simmetrica
if isequal(A, A')
    disp('A simmetrica');
    nA = size(A, 1);

    % la definita positività significa che x'Ax >= 0 per ogni x
    % appartenente ad R^n, di fatto non verificabile con la definizione:
    % applichiamo il criterio di Sylvester, secondo cui affinché A sia DP
    % il determinante di tutte le sottomatrici deve essere positivo

    % calcolo il determinante delle sottomatrici (criterio di Sylvester)
    for ii = 1:nA
        if det(A(1:ii, 1:ii)) > 0
            % se siamo arrivati all'ultima iterazione allora le condizioni
            % su tutti i determinanti sono verificate
            if ii == nA
                disp('A DP');
            end
        else
            error('A non DP');
        end
    end
else
    disp('A non simmetrica');
end

%% 1.2 B a dominanza diagonale stretta?

% estraggo la diagonale principale e la dimensione
dB = diag(B);
nB = size(B, 1);

% per colonne
% inizializzo un vettore riga vuoto
r = zeros(1, nB);
for jj = 1:nB
    % sommo gli elementi extradiagonali (con indice di riga diverso da jj)
    % sulle varie colonne e salvo in un vettore
    r(jj) = sum(abs(B([1:jj-1, jj+1:nB], jj)));
end
% if su un vettore è vero solo se le condizioni son verificate per tutti i
% suoi elementi
if dB - r' > 0
    disp('B a dominanza diagonale per colonne');
else
    disp('B non a dominanza diagonale per colonne');
end

% per righe
% inizializzo un vettore colonna vuoto
r = zeros(nB, 1);
for ii = 1:nB
    % sommo gli elementi extradiagonali (con indice di colonna diverso da
    % jj) sulle varie righe e salvo in un vettore
    r(ii) = sum(abs(B(ii, [1:ii-1, ii+1:nB])));
end
% if su un vettore è vero solo se sono verificate tutte le condizioni
if dB - r > 0
    disp('B a dominanza diagonale per righe');
else
    disp('B non a dominanza diagonale per colonne');
end

%% 1.3 C condizione necessaria e sufficiente?

% inizializzo l'indice, la dimensione e determinante arbitriario non nullo
ii = 1;
nC = size(C, 1);
dC = 1;

% finché non sono all'ultimo minore e il determinante del precedente non è
% nullo (si potrebbe fare con un ciclo for)
while ii < nC && dC ~= 0
    % verifico la condizione e incremento l'indice
    dC = det(C(1:ii, 1:ii));
    ii = ii + 1;
end

% utlizzando un ciclo for
% for ii = 1:nC-1
%     dC = det(C(1:ii, 1:ii));
%     if dC == 0
%         break
%     end
% end

% se sono uscito dal ciclo perché ho finito è LU-fattorizzabile
if dC ~= 0
    disp('C LU-fattorizzabile')
else
    disp('C non LU-fattorizzabile')
end

%% 2. calcolo fattorizzazione LU

% per il nostro corso gli facciamo sempre fare il pivoting, chiediamo
% sempre la marice P negli argument della funzione lu()

[LA, UA, PA] = lu(A);

% verifichiamo che effettivamente siano uguali
norm(PA*A - LA*UA)

% per verificare se è stato effettuato pivoting
if isequal(PA, eye(size(A)))
    disp("Pivoting non effettuato");
else
    disp("Pivoting effettuato");
end

[LB, UB, PB] = lu(B);

[LC, UC, PC] = lu(C);

%% 3. sostituzione in avanti e all'indietro

% creiamo due funzioni che permettono di effettuare le sostituzioni in
% avanti e all'indietro
% vedi file fwsub.m e bksub.m

%% 4. soluzione Ax = b

% il testo ci assegna la soluzione esatta x_ex da cui ricavarci il temine
% noto b
x_ex = ones(nA, 1);
b = A * x_ex;

% risolviamo il sistema lineare e troviamo la soluzione approssimata al
% calcolatore x
y = fwsub(LA, PA * b);
x = bksub(UA, y);

%% 5. calcolo della norma 2 dell'errore relativo

% calcoliamo l'errore relativo e il residuo relativo
err_rel = norm(x - x_ex) / norm(x_ex);
res_norm = norm(b - A*x)/norm(b);

% calcoliamo in numero di condizionamento della matrice A
cond(A);
% questa matrice ha un numero di condizionamento basso, è ragionevole
% ritrovare un errore ridotto