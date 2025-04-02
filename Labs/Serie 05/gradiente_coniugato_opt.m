function [xvect, it] = gradiente_coniugato_opt(Phi, GradPhi, flag_metodo, x0, toll, nmax )

% function [xvect, it] = gradiente_coniugato(Phi, GradPhi, flag_metodo, x0, toll, nmax )
% algoritmo gradiente (gradiente coniugato) per la risoluzione di problemi 
% di ottimizzazione per funzionali in 2-dimensioni. Tramite flag_metodo è
% possibile scegliere se usare il metodo del gradiente o il metodo del
% gradiente coniugato (con \beta_k di Fletcher-Reeves)
% IN
%   - Phi: funzionale inline/anonymus
%   - GradPhi: gradiente di Phi come funzione inline/anonymus
%   - flag_metodo: flag per scegliere tra gradiente e gradiente coniugato
%       - 'G': gradiente
%       - 'FR': gradiente coniugato (Fletcher-Reeves)
%   - x0: guess iniziale
%   - toll: tolleranza
%   - nmax: numero massimo di iterazioni
% OUT
%   - xvect: matrice che su ogni colonna contiene la soluzione a ogni
%       iterazione (l'ultima colonna è la soluzione)
%   - it: numero di iterazioni

% inizializzo variabili per ciclo iterativo
it = 0;
err = toll+1;

% inizializzo la matrice soluzione
xvect = x0;

% calcolo la direzione iniziale di discesa
d = -GradPhi(x0(1), x0(2));

% parametri backtracking
c1_bt = 1e-4;
rho_bt = 0.3;
nmax_bt = 10;

% ciclo while per calcolo di x^{it+1} (it = 0, 1, ...)
while it < nmax && err > toll
    
    % calcolo il nuovo valore di alpha_k
    [alpha_k, ~] = backtracking(Phi, GradPhi, x0, d, c1_bt, rho_bt, nmax_bt);
     
    % calcolo il nuovo valore di x
    x_new = x0 + alpha_k * d;

    % swtich per scelta metodo gradiente o gradiente coniugato
    switch flag_metodo
        case 'G' % gradiente
            beta_k = 0;
        case 'FR' % gradiente coniugato (Fletcher-Reeves)
            beta_k = norm(GradPhi(x_new(1), x_new(2)))^2 / norm(GradPhi(x0(1), x0(2)))^2;
        otherwise
            error('metodo CG non definito');
    end

    % aggiorno la direzione di discesa (se beta_k = 0 -> gradiente)
    d = -GradPhi(x_new(1), x_new(2)) + beta_k * d;

    % aggiorno quantità per metodo iterativo (it, err, x0, xvect)
    it = it + 1;
    err = norm(GradPhi(x0(1), x0(2)));
    x0 = x_new;
    xvect = [xvect, x_new];
end

% printo informazioni sulla convergenza o meno
switch flag_metodo
    case 'G' % gradiente
        if err <= toll
             fprintf(['\n Il metodo del gradiente converge in %d iterazioni \n' ...
                 'al punto di minimo x = (%f,%f)^T \n '], it, xvect(1,end), xvect(2,end));
        else
             fprintf('\n Il metodo del gradiente non converge in %d iterazioni. \n', it)
        end
        
    case 'FR' % gradiente coniugato (Fletcher-Reeves)
        if err <= toll
             fprintf(['\n Il metodo del gradiente coniugato (Fletcher-Reeves) converge in %d iterazioni \n' ...
                 'al punto di minimo x = (%f,%f)^T \n '], it, xvect(1,end), xvect(2,end));
        else
             fprintf('\n Il metodo del gradiente coniugato (Fletcher-Reeves) non converge in %d iterazioni. \n',  it)
        end    
end

end