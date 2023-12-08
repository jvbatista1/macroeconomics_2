% Parâmetros
beta = 0.96;
alpha = 0.25;
A = (alpha * beta)^(-1);

% Parâmetros para k
n = 41;
kappa = 3.5 / (n - 1);
k_values = kappa * (0:n-1) + kappa;

% Parâmetros para z
sigma = 0.5;
m = 41;
delta_z = (6 * sigma) / (m - 1);
z_values = -3 * sigma + (0:m-1) * delta_z;

% Calcula q_j
phi = @(x) exp(-0.5 * x.^2) / sqrt(2 * pi);
q_values = phi(z_values / sigma) / sum(phi(z_values / sigma));
q = q_values';

% Inicialização
V = zeros(n, m);
epsilon = 1e-5;
max_diff = inf;

% Iteração da função valor
while max_diff > epsilon
    % Passo 2.1: Crie a matriz de Utilidade U
    U = zeros(n * m, n);
    for j = 1:m
        for i = 1:n
            U((i - 1) * m + j, :) = log(exp(z_values(j)) * A * k_values(i)^alpha - k_values);
        end
    end
    
    % Passo 2.2: Compute EV^l
    %EV = reshape(V * q, n, m);
    %EV = reshape(V(:) * q, n, m);
    %EV = reshape(V .* q', n, m);
    %EV = V * q;
    %EV = reshape(EV, n, m);
    %EV = reshape(V * q, m, n)';
    %EV = reshape(sum(V .* q, 2), n, m);
    EV = V * q;
    
    % Passo 2.3: Compute W^{l+1}
    W = U + beta * kron(ones(n, 1), EV);

    % Passo 2.4: Compute V^{l+1}
    V_new = reshape(max(W, [], 2), n, m);
    
    % Verifique a condição de parada
    max_diff = max(abs(V_new(:) - V(:)));
    
    % Atualize V
    V = V_new;
end

% Encontrar a função política
P = zeros(n * m, m);
for j = 1:m
    [~, max_index] = max(U + beta * kron(ones(m, 1), EV(:)), [], 2);
    P((1:n) + (j - 1) * n, max_index) = 1;
end

% Exibir resultados
disp('Função Valor Convergida:');
disp(V);
disp('Função Política:');
disp(P);
