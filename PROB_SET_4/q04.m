%%%%%%%%%%%%
% Macroeconomia II
% Lista 4
% Questão 4 (Huggett 1993)
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%

% Parâmetros
beta = 0.96;
gamma = 2;
w = 1;
phi = -2;

% Função Valor Iterativa
[V, A] = iterate_value_function(beta, gamma, w, phi);

% Método da Bisseção para encontrar a taxa de juros de equilíbrio
r_low = 0.01;
r_high = 0.10;
r_equilibrium = find_equilibrium_interest_rate(beta, gamma, w, phi, r_low, r_high);

% Simulação da dinâmica da economia
num_agents = [100, 1000, 10000];
num_periods = [100, 1000, 10000];

for i = 1:length(num_agents)
    for j = 1:length(num_periods)
        num_simulations = num_agents(i);
        num_time_periods = num_periods(j);
        
        % Simulação da dinâmica da economia
        [asset_distribution, total_assets_over_time] = simulate_economy(beta, gamma, w, phi, r_equilibrium, num_simulations, num_time_periods);
        
        % Plote o estoque total de ativos
        figure;
        plot(1:num_time_periods, total_assets_over_time);
        title(['Dinâmica da Economia - ', num2str(num_simulations), ' Agentes, ', num2str(num_time_periods), ' Períodos']);
        xlabel('Período de Tempo');
        ylabel('Estoque Total de Ativos');
    end
end

% Histograma da distribuição de ativos no último período
final_asset_distribution = asset_distribution(:, end);
figure;
hist(final_asset_distribution, 50);
title('Histograma da Distribuição de Ativos no Último Período');
xlabel('Ativos');
ylabel('Frequência');
