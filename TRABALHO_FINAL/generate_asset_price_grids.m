function [k_grid, p_grid] = generate_asset_price_grids()

    % Parâmetros adicionais
    beta_val = 0.96; % Substitua pelo valor real

    % Calcula Kss
    alpha = 0.36;
    delta = 0.08;
    rss = 1/beta_val - 1;
    Kss = (rss + delta) / alpha ^ (1 / (alpha - 1));

    % Número de pontos no grid de ativos
    n_k = 512;

    % Criar grid de ativos
    k_grid = zeros(1, n_k);
    % Um terço dos pontos no intervalo [0, Kss]
    k_grid(1:floor(n_k/3)) = linspace(0, Kss, floor(n_k/3));
    % Um terço dos pontos no intervalo (Kss, 3Kss]
    k_grid(floor(n_k/3)+1:2*floor(n_k/3)) = linspace(Kss, 3*Kss, floor(n_k/3));
    % Um terço dos pontos no intervalo (3Kss, 15Kss]
    k_grid(2*floor(n_k/3)+1:end) = linspace(3*Kss, 15*Kss, n_k-2*floor(n_k/3));

    % Número de pontos no grid de preços
    n_p = 251;

    % Criar grid de preços
    p_grid = zeros(1, n_p);
    % Um terço dos pontos no intervalo [-delta, 0)
    p_grid(1:floor(n_p/3)) = linspace(-delta, 0, floor(n_p/3));
    % Dois terços dos pontos no intervalo [0, rss]
    p_grid(floor(n_p/3)+1:end) = linspace(0, rss, n_p-floor(n_p/3));

end

