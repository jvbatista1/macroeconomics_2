function simulation_results = simulate_model()

    % Parâmetros do Modelo
    beta = 0.96;
    beta_val = 0.96;
    alpha = 0.36;
    delta = 0.08;
    mu_values = [1, 3, 5];
    rho_values = [0, 0.3, 0.6, 0.9];
    sigma_values = [0.2, 0.4];
    
    % Configuração de Simulação
    k_grid_size = 512;
    h_values = linspace(0.8, 1.2, 10);
    p_grid_size = 251;

    % Iteração da Função Valor
    max_iter = 1000;
    tolerance = 1e-5;

    % Inicialização da Função Valor e Política de Poupança
    V = zeros(k_grid_size, length(h_values), length(mu_values), length(rho_values), length(sigma_values), length(p_grid_size));
    k_prime_policy = zeros(size(V));

    % Outros Parâmetros
    beta = 0.96;
    tau = 0.1; % Substitua pelo valor real
    [z_values, transition_matrix] = tauchen_method(0, 0.2, 21, 3);
    [k_grid, p_grid] = generate_asset_price_grids();

    % Implementação da iteração da função valor
for mu_index = 1:length(mu_values)
    for rho_index = 1:length(rho_values)
        for sigma_index = 1:length(sigma_values)
            mu = mu_values(mu_index);
            rho = rho_values(rho_index);
            sigma = sigma_values(sigma_index);

           % Inicialização da Função Valor
           V(:, :, mu_index, rho_index, sigma_index, :) = repmat((log(k_grid) - tau * k_grid) / (1 - beta), length(h_values), length(p_grid));

            % Iteração da Função Valor
            for iter = 1:max_iter
                V_old = V(:, :, mu_index, rho_index, sigma_index);

                % Atualização da Função Valor
                for h_index = 1:length(h_values)
                    h = h_values(h_index);
                    for k_index = 1:k_grid_size
                        k = k_grid(k_index);
                        for p_index = 1:p_grid_size
                            p = p_grid(p_index);

                            % Restrição de recursos
                            c = max(w * h + (1 + r - p) * k - k_grid, 0);

                            % Utilidade
                            u = log(c) - tau * c;

                            % Máximo da equação de Bellman
                            [max_value, max_index] = max(u + beta * sum(sum(V(:, :, mu_index, rho_index, sigma_index, :), 2) .* transition_matrix(h_index, :)));

                            % Atualização da Função Valor
                            V_new(k_index, h_index, mu_index, rho_index, sigma_index, p_index) = max_value;
                        end
                    end
                end

                % Verificar Convergência
                if max(abs(V_new(:) - V_old(:))) < tolerance
                    break;
                end

                V(:, :, mu_index, rho_index, sigma_index) = V_new;
            end
        end
    end
end


    % Agora, V contém os valores da função valor após a iteração

    % Você pode prosseguir com a análise, por exemplo, derivar a política de poupança
    % k_prime = g(k, h), calcular os agregados e assim por diante.

    % Armazenar resultados em uma estrutura
    simulation_results.V = V;
    simulation_results.k_prime_policy = k_prime_policy;

end