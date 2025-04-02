function plot_phi(x, y, Phi, title)

% plot_phi(x, y, Phi, title)
% Plot della funzione Phi (funzione in due variabili) con comandi surf e  
% contourf
% 
% INPUT x     : (double, vettore) coordinate asse x
%       y     : (double, vettore) coordinate asse y
%       Phi   : (handle function) Funzionale
%       title : (string) titolo del grafico (opzionale)

if nargin == 4
    sgtitle(title);
end

[X,Y] = meshgrid(x,y); 
PHI   = Phi(X,Y);

subplot(1,2,1)
surf(X,Y,PHI,'Lines','no')
hold on;
xlabel('x');
ylabel('y');
zlabel('z')

subplot(1,2,2)
contourf(X, Y, PHI, 30, 'LineWidth', 0.5,'FaceAlpha',1);
hold on;
axis equal;
colorbar;
xlabel('x');
ylabel('y');

end

