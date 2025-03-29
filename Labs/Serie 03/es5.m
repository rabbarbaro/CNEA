clc
clear
close all

%% 1

% importiamo l'immagine 'street1.jpg' implementata in matlab e la
% trasformiamo in scala di grigi
IM = imread('street1.jpg');
IM = rgb2gray(IM);
imshow(IM)
% l'immagine è in interi, convertiamo (facciamo un cast) in virgola mobile
A = double(IM);

[m, n] = size(A);

%% 2

% calcoliamo la SVD di A
[U, S, V] = svd(A);

% verifichiamo che è effettivamente la decomposizione di A
% norm(U * S * V' - A)
% la norma dell'errore è circa zero

p = min(n, m);
% il rango della matrice (in generale) è il minimo delle due dimensioni

% scegliamo k come da testo e assembliamo Ak, moltiplicando i termini della
% SVD ristretti al rango k scelto
k = 10;
Ak = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';

%% punto 3

close all

% plottiamo l'immagine non compressa e quella compressa

figure
subplot(2,2,1)
imshow(IM)
title('Uncompressed image')

subplot(2,2,3)
imshow(uint8(Ak))
title('Compressed image')

% in scala logaritmica plotto il valore dei valori singolari della matrice
% rispetto alla loro posizione (il loro contributo)
% aggiungo una linea tratteggiata dove tronchiamo
subplot(1,2,2)
loglog(1:p, diag(S), 'b-', [k k], [S(p,p) S(1,1)],'r--','LineWidth',2)
grid on
xlabel('i')
ylabel('\sigma_i')
title('\sigma_i vs. i')

%% punto 4

% calcoliamo le dimensioni non copressa e compressa, e poi il rapporto di
% compressione

uncompressed_size = m*n;
compressed_size = k * (m+n+1);
compression_ratio = uncompressed_size/compressed_size;

fprintf("Uncompressed size: %d\n", uncompressed_size)
fprintf("Compressed size: %d\n", compressed_size)
fprintf("Compression ratio: %f\n", compression_ratio)

% calcoliamo l'errore relativo sapendo dalla teoria (TH di Schmidt-Eckart-
% Young) come varia l'errore relativo in norma di Frobenius
relative_error = sqrt(sum(diag(S(k+1:end, k+1:end)).^2)) / ...
    sqrt(sum(diag(S).^2));
fprintf("Relative error: %f\n", relative_error)