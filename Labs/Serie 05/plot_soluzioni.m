function plot_soluzioni(Phi, xvect, range_xlim, range_ylim, x_min_ex)

% plot_soluzioni(Phi, xvect, range_xlim, range_ylim, x_min_ex)
% Plot delle successioni di soluzione algoritmo di ottimizzazione. Viene
% consigliato l'utilizzo dopo aver eseguito plot_phi(x, y, Phi, title).
% 
% INPUT Phi        : (handle function) Funzionale
%       xvect      : (double, matrice) matrice che, su ogni colonna, contiene 
%                    la soluzione a ogni iterazione 
%       range_xlim : (double, vettore) limiti asse x
%       range_ylim : (double, vettore) limiti asse y
%       x_min_ex   : (double, vettore) soluzione esatta (opzionale)

subplot(1, 2, 1);
plot3( xvect(1,:),xvect(2,:), Phi(xvect(1,:),xvect(2,:)),'ko--','MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
if(nargin == 5)
    plot3(x_min_ex(1), x_min_ex(2), Phi(x_min_ex(1), x_min_ex(2)), "o",'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 1.5);
end
subplot(1, 2, 2);
plot( xvect(1,:),xvect(2,:),'ko--','MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
if(nargin == 5)
    plot(x_min_ex(1), x_min_ex(2), "o",'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'LineWidth', 1.5);
end
colorbar;
axis equal;
xlabel('x');
ylabel('y');
xlim(range_xlim);
ylim(range_ylim);

end

