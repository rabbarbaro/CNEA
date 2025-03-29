clear
clc
close all

% inizializzo l'anno 0, definico il deposito inizale e il tasso
anni = 0;
euri = 1e4;
interesse = 1.02;

% finché sono sotto al milione aggiungo un interesse del
while euri < 1e6
    euri = interesse * euri + 1e4;
    anni = anni + 1;
end