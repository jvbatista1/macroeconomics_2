%%%%%%%%%%%%
% Macroeconomia II
% Lista 4
% Questão 2 (Huggett 1993)
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

% Número de simulações
N = 1000;

% Inicialização
s_history = zeros(N, 1);
s_history(1) = S(randi(length(S))); % Escolher um estado inicial aleatório

% Simulação da cadeia de Markov
for t = 2:N
    % Gerar número aleatório entre 0 e 1
    random_number = rand();
    
    % Encontrar o novo estado com base nas probabilidades cumulativas
    cumulative_probabilities = cumsum(P(s_history(t-1)==S, :));
    new_state_index = find(random_number <= cumulative_probabilities, 1, 'first');
    s_history(t) = S(new_state_index);
end

% Calcular persistência
persistence = sum(s_history(1:end-1) == s_history(2:end)) / (N - 1);

% Calcular variância
variance = var(s_history);

% Exibir resultados
disp('Histórico das realizações da cadeia de Markov:');
disp(s_history);
fprintf('Persistência: %.4f\n', persistence);
fprintf('Variância no longo prazo: %.4f\n', variance);
