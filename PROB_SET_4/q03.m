%%%%%%%%%%%%
% Macroeconomia II
% Lista 4
% Questão 3 (Huggett 1993)
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%

% Definir a matriz de transição
P = [0.7497, 0.2161, 0.0322, 0.002, 0;
     0.2161, 0.4708, 0.2569, 0.0542, 0.002;
     0.0322, 0.2569, 0.4218, 0.2569, 0.0322;
     0.002, 0.0542, 0.2569, 0.4708, 0.2161;
     0, 0.002, 0.0322, 0.2161, 0.7497];

% Vetor de estados
S = [0.6177, 0.8327, 1.0000, 1.2009, 1.6188];

% Inicialização
epsilon = 1e-5; % Critério de parada
xi_prime = ones(length(S), 1) / length(S); % Inicia com distribuição uniforme

% Iteração para encontrar a distribuição estacionária
while true
    xi = xi_prime;
    xi_prime = xi' * P;
    
    % Verificar o critério de parada
    if norm(xi_prime - xi) < epsilon
        break;
    end
end

% Calcular oferta total de trabalho
N = sum(xi .* S);

% Exibir resultados
disp('Distribuição Estacionária:');
disp(xi');
fprintf('Oferta Total de Trabalho: %.4f\n', N);
