%%%%%%%%%%%%
% Macroeconomia II
% Lista 2
% Questão 2
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%

% Parâmetros
beta = 0.9;
c = 0.2;
salarios = 1:100;
probabilidades = ones(1, 100) * 0.01;

% Inicialização da função valor
V = zeros(size(salarios));

% Critério de convergência
tolerancia = 1e-5;
diferenca = inf;

% Iteração de Bellman
while diferenca > tolerancia
    V_anterior = V;

    % Calcular novo valor para cada salário
    for i = 1:length(salarios)
        salario_atual = salarios(i);
        
        % Calcular valor esperado futuramente
        valor_esperado_futuro = sum(probabilidades .* V);
        
        % Atualizar função valor usando a equação de Bellman
        V(i) = max(salario_atual / (1 - beta), c + beta * valor_esperado_futuro);
    end

    % Calcular diferença entre iterações consecutivas
    diferenca = max(abs(V - V_anterior));
end

% Exibir resultado
disp('Função Valor:')
disp(V);

%% item 2
% Parâmetros
beta = 0.9;
c = 0.2;
salarios = 1:100;

% Calcular a função de densidade de probabilidade
probabilidades = 1 ./ salarios;
probabilidades = probabilidades / sum(probabilidades);

% Inicialização da função valor e do salário de reserva
V = zeros(size(salarios));
salario_reserva = zeros(size(salarios));

% Critério de convergência
tolerancia = 1e-5;
diferenca = inf;

% Iteração de Bellman
while diferenca > tolerancia
    V_anterior = V;

    % Calcular novo valor para cada salário
    for i = 1:length(salarios)
        salario_atual = salarios(i);
        
        % Calcular valor esperado futuramente
        valor_esperado_futuro = sum(probabilidades .* V);
        
        % Atualizar função valor usando a equação de Bellman
        V(i) = max(salario_atual / (1 - beta), c + beta * valor_esperado_futuro);
        
        % Verificar se o trabalhador é indiferente
        if abs(salario_atual / (1 - beta) - (c + beta * valor_esperado_futuro)) < tolerancia
            salario_reserva(i) = salario_atual;
        end
    end

    % Calcular diferença entre iterações consecutivas
    diferenca = max(abs(V - V_anterior));
end

% Exibir resultado
disp('Função Valor:')
disp(V);

disp('Salário de Reserva:')
disp(salario_reserva);

