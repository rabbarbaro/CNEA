clear
clc
close all

%% 1

% costruisco A1 come da testo
A1 = [6.8 -2.4
    -2.4 8.2];

% verifico che abbia effettivamente gli autovalori e autovettori indicati
[v, l] = eig(A1);
fprintf("A ha autovalori:")
diag(l)
fprintf("A ha aurovettori:")
v

% verifico a mano autovettori e autovalori
% v1 = 0.2 * [4 3]'
% A1 * v1
% l1 = norm(A1 * v1) / norm(v1)
% v2 = 0.2 * [-3 4]'
% A1 * v2
% l2 = norm(A1 * v2) / norm(v2)

% posso costruire la matrice A1 anche a partire da autovettori e autovalori
lambda(1) = 5;
lambda(2) = 10;
v1 = 0.2 * [4 3]';
v2 = 0.2 * [-3 4]';
V  = [v1,v2]; 
D1 = diag(lambda);
A1 = V * D1 * V';

% genero anche A2 in questo modo
lambda(1) = 1;
D2 = diag(lambda);
A2 = V * D2 * V';

% definisco b come dal testo
b = [4 8]';

% costruisco la mesh per fare il grafico 3D
[X,Y] = meshgrid(-10:0.5:10,-10:0.5:10);

% definisco le due forme quadratiche dalle matrici assegnate
phi1 = @(X,Y) 0.5 * (A1(1,1)*X.^2 + 2*A1(1,2)*X.*Y + A1(2,2)*Y.^2) - (b(1)*X + b(2)*Y);
phi2 = @(X,Y) 0.5 * (A2(1,1)*X.^2 + 2*A2(1,2)*X.*Y + A2(2,2)*Y.^2) - (b(1)*X + b(2)*Y);

% plotto la superficie e le linee di livello delle due forme quadratiche
subplot(2,2,1)
surf(X,Y,phi1(X, Y),'Lines','no');
title('\phi_1')
subplot(2,2,2)
contour(X,Y,phi1(X, Y), 15)
title('Linee di livello \phi_1')
subplot(2,2,3)
surf(X,Y,phi2(X, Y),'Lines','no');
title('\phi_2')
subplot(2,2,4)
contour(X,Y,phi2(X, Y), 15)
title('Linee di livello \phi_2')

%% 2

% vedi file richardson_it.m

% definisco la guess iniziale, il numero di iterazioni massime e la
% tolleranza
x0 = [-9 -9]';
nmax = 1000;
toll = 1e-6;

% alpha1 come da testo
alpha1 = 0.05;

% chiamo richardson_it mettendo P = identità
[x_rich, n_it_rich] = richardson_it(A2, b, eye(2), x0, toll, nmax, alpha1);

% plottiamo
figure
subplot(1, 3, 1)
contour(X,Y,phi2(X, Y), 50)
title('\alpha = 0.05')
axis equal
hold on
plot(x_rich(1,:),x_rich(2,:),'-or')

% alpha2 come da testo
alpha2 = 0.24;

% chiamo richardson_it mettendo P = identità
[x_rich2, n_it_rich2] = richardson_it(A2, b, eye(2), x0, toll, nmax, alpha2);

% plottiamo
subplot(1, 3, 2)
contour(X,Y,phi2(X, Y), 50)
title('\alpha = 0.24')
axis equal
hold on
% dato che esplode limito le iterazioni a 10
plot(x_rich2(1,1:10),x_rich2(2,1:10),'-or')

% per il gradiente basta non imporre alpha, sempre P = I
[x_grad,n_it_grad] = richardson_it(A2, b, eye(2), x0, toll, nmax);

% plottiamo
subplot(1, 3, 3)
contour(X,Y,phi2(X, Y), 50)
title('Gradiente')
axis equal
hold on
plot(x_grad(1,:),x_grad(2,:),'-or')

%% 3

% definisco il precondizionatore
P = [1.0912 -0.8587
    -0.8587 1.5921];

% calcolo P^-1A e P^-1b
Aprec = P \ A2;
bprec = P \ b;

% definisco la nuova forma quadratica
phiprec = @(X,Y) 0.5 * (Aprec(1,1)*X.^2 + 2*Aprec(1,2)*X.*Y + Aprec(2,2)*Y.^2) - (bprec(1)*X + bprec(2)*Y);

% plotto le due forme quadratiche
figure
subplot(1,2,1)
contour(X,Y,phi2(X, Y), 15)
title('Linee di livello \phi_2')
subplot(1,2,2)
contour(X,Y,phiprec(X, Y), 15)
title('Linee di livello \phi_prec')
% ha lo stesso minimo ma le curve di livello sono meno schiacciate

% chiamiamo il metodo del gradiente non specificando alpha (utilizziamo
% il metodo precondizionato con P)
[x_prec,n_it_prec] = richardson_it(A2, b, P, x0, toll, nmax);

% plottiamo
figure;
contour(X,Y,phi2(X, Y), 50)
axis equal
hold on
title('Gradiente Precondiz. vs Gradiente')
% in grigio le direzioni del gradiente normale
plot(x_grad(1,:),x_grad(2,:),'--','Color',[0.48 0.48 0.48])
% in rosso quelle del gradiente precondizionato
plot(x_prec(1,:),x_prec(2,:),'-or')

%% 4

% chiamiamo il metodo del gradiente coniugato non precondizionato
[x_cg, ~, ~, n_it_cg] = pcg(A2, b, toll, nmax, eye(size(A2)), eye(size(A2)), x0);
fprintf('Il metodo del gradiente coniugato converge in %d iterazioni', n_it_cg);

% effettivamente dalla teoria sappiamo che converge alla soluzione esatta
% al massimo in n iterazioni