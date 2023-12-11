function [V, A] = iterate_value_function(beta, gamma, w, phi)
    % Parâmetros do modelo
    num_states = 5; % Número de estados no vetor S
    epsilon = 1e-6; % Critério de convergência

    % Inicialização
    V = zeros(num_states, 1);
    A = zeros(num_states, 1);
    delta = epsilon + 1;

    while delta > epsilon
        V_new = zeros(num_states, 1);
        A_new = zeros(num_states, 1);

        for i = 1:num_states
            % Código para calcular V_new(i) e A_new(i) usando a função utilidade
            % e a equação de Bellman.

            % Exemplo (pode não ser preciso para o seu modelo):
            % V_new(i) = max_utility(i); 
            % A_new(i) = argmax_utility(i);
        end

        % Atualização
        delta = max(abs(V_new - V));
        V = V_new;
        A = A_new;
    end
end