clear
clc
close all

% Carica l'immagine e convertila in scala di grigi
img = imread('cat.png');  % Sostituisci con il nome della tua immagine
img_gray = rgb2gray(img);  % Converti in scala di grigi se necessario
A = double(img_gray);  % Converti in formato double per la SVD

% Calcola la decomposizione SVD
[U, S, V] = svd(A);

% Valori di k per la ricostruzione
k_values = [5, 20, 50, 100];  

[m,n] = size(A);
uncompressed_size = m*n;

% Mostra le immagini ricostruite
figure;
subplot(2,3,1);
imshow(uint8(A));  
title('Immagine originale');
xlabel(["Uncompressed size: " num2str(uncompressed_size)]);

for i = 1:length(k_values)
    k = k_values(i);
    
    % Approssimazione con i primi k valori singolari
    A_k = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';

    compressed_size = k * (m+n+1);
    
    compression_ratio = uncompressed_size/compressed_size;
    relative_error = sqrt(sum(diag(S(k+1:end,k+1:end)).^2)) / ...
            sqrt(sum(diag(S).^2));

    % Mostra l'immagine compressa
    subplot(2,3,i+1);
    imshow(uint8(A_k));  
    title(['k = ' num2str(k)]);
    xlabel(["Compression ratio: " num2str(compression_ratio) ...
        "Relative error: " num2str(relative_error)]);
end
