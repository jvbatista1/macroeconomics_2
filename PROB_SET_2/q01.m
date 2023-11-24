%%%%%%%%%%%%
% Macroeconomia II
% Lista 2
% Questão 1
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%

%% item ii
% Parâmetros
beta = 0.96;
alpha = 0.25;
kappa = 0.01;
n = 41;
k1 = 0.8;
kn = 1.2;
A = (alpha * beta)^(-1);

% Matriz de utilidade U
k = linspace(k1, kn, n);
U = zeros(n, n);

for i = 1:n
    for j = 1:n
        U(i, j) = log(A * k(i)^alpha - k(j));
    end
end

% Vetor tildeV inicial (pode ser inicializado com zeros)
tildeV = zeros(n, 1);

% Iterações
num_iterations = 1000;

for l = 1:num_iterations
    % Computar matriz W^(l+1)
    ones_vector = ones(n, 1);
    tildeV_transposed = tildeV';
    W = U + beta * kron(ones_vector, tildeV_transposed);
    
    % Atualizar tildeV^(l+1)
    tildeV_new = max(W, [], 2);
    
    % Condição de convergência
    if max(abs(tildeV_new - tildeV)) < 1e-6
        break;
    end
    
    tildeV = tildeV_new;
end

% Plotagem dos resultados
figure;
plot(k, tildeV, '-o');
xlabel('k');
ylabel('\tilde{V}');
title('Valor Aproximado \tilde{V}');

%% item iii
% Vetor \tilde{V}*
tildeV_star = tildeV;

% Computar matriz de política P
P = zeros(n, n);

for i = 1:n
    [~, j_star] = max(U(i, :) + beta * tildeV_star');
    P(i, j_star) = 1;
end

% Teste da política multiplicando pelo vetor de capital atual
k_j_star = P * k';

%% item 2
% Plotagem dos resultados
figure;
plot(k, k_j_star, '-o');
xlabel('k_i');
ylabel('k_{j^*}');
title('Função de Política Ótima');

% Função de política real
k_real_policy = alpha * beta * A * k.^alpha;

% Plotagem dos resultados
figure;
hold on;

% Gráfico da função de política obtida pela iteração de valor
plot(k, k_j_star, '-o', 'DisplayName', 'Política Iterativa');

% Gráfico da função de política real
plot(k, k_real_policy, '-+', 'DisplayName', 'Política Real');

xlabel('k_i');
ylabel('k_{j^*}');
title('Comparação da Função de Política');
legend('show');

hold off;