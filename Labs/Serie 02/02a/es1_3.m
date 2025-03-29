clear
clc
close all

% abbiamo 2(n+1) incognite: V_k e I_k per k = 1:n, I_0, V_0
% - I0 è assegnata (-1)
% - essendo le n resistenze e il generatore di corrente in parallelo la V
%   ai loro capi sarà la stessa (-n)
% - per la OL ogni resistenza avrà V_k = V = R_k * I_k (-1, toglie la
%   dipendenza da V che non è richiesta)
% otteniamo alla fine n incognite
% abbiamo quindi bisogno di n equazioni:
% - per la KCL I0 = I_1 + I_2 + ... + I_n (1)
% - per la KVL ogni maglia avrà R_k * I_k = R_k-1 * I_k-1 per k = 2:n (n-1)

% il sistema lineare risultante sarà descritto da (ogni elemento su ogni
% colonna va moltiplicato per l'incognita corrispontente I_1 ... I_n):
% A = [1    1     1  ...   1     1
%     R_1 -R_2    0  ...   0     0
%      0   R_1  -R_3 ...   0     0
%     ...
%      0    0     0  ... R_n-1 -R_n]
% b = [I_0 0 0 ... 0]'

%% 1

% definisco i dati come dal testo
n = 20;
R = ones(n,1);
I0 = 2;
% definisco b come dal testo
b = zeros(n, 1);
b(1) = I0;
% definisco A come dal testo
A = -diag(R) + diag(R(1:n-1), -1);
A(1, :) = 1;

%% 2

% calcolo la fattorizzazione LU di A con pivoting per righe
[LA, UA, PA] = lu(A);

% verifico che non è stato effettuato pivoting
if isequal(PA, eye(size(A)))
    disp("Pivoting non effettuato");
else
    disp("Pivoting effettuato");
end

%% 3

% risolvo numericamente con fwsub.m e bksub.m il sistema Ax = b
y = fwsub(LA, PA*b);
x = bksub(UA, y);

%% 4

% calcoliamo norma 2 dell'errore relativo e del residuo relativo avendo
% assegnta dal testo la soluzione esatta x_ex
x_ex = (I0 / n) * ones(n, 1);
err_rel = norm(x - x_ex) / norm(x_ex);
res_norm = norm(b - A*x)/norm(b);

% calcolo il numero di condizionamento di A
cond(A);
% A ha un numero di condizionamento basso, è ragionevole ritrovare un
% errore ridotto

%% 5

% modifico R_1 come dal testo (R_1 compare solo in una posizione in A)
A(2,1) = 1e3;

% calcolo nuovamente la soluzione numerica
[LA, UA, PA] = lu(A);
y = fwsub(LA, PA*b);
x = bksub(UA, y);

% verifico che questa volta è stato effettuato pivoting
if isequal(PA, eye(size(A)))
    disp("Pivoting non effettuato");
else
    disp("Pivoting effettuato");
end

% calcolo il nuovo residuo relativo
res_norm = norm(b - A*x)/norm(b);
% calcolo il nuovo numero di condizionamento di A
cond(A);
% A ha un numero di condizionamento relativamente alto
% il residuo rimane comunque basso, ma l'errore sarà più elevato