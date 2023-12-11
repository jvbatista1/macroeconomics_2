%%%%%%%%%%%%
% Macroeconomia II
% Lista 4
% Questão 1 (Huggett 1993)
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%

% Definir a matriz de transição
P = [0.5, 0.5, 0;
     0.5, 0.2, 0.3;
     0, 0.7, 0.3];

% Encontrar a distribuição estacionária
[V, D] = eig(P');  % Calcular autovalores e autovetores transpostos
[~, idx] = max(abs(diag(D)));  % Encontrar índice do autovalor dominante
pi_estacionaria = abs(V(:, idx) / sum(V(:, idx)));  % Normalizar o autovetor associado

disp('Distribuição Estacionária:');
disp(pi_estacionaria);
