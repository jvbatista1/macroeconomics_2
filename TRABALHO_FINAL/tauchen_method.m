function [z_values, transition_matrix] = tauchen_method(rho, sigma, n_z, m)

% Parâmetros do Método de Tauchen
z_values = linspace(-m * sqrt(sigma^2 / (1 - rho^2)), m * sqrt(sigma^2 / (1 - rho^2)), n_z);

% Intervalo entre os estados
z_interval = z_values(2) - z_values(1);

% Probabilidade de transição
transition_matrix = zeros(n_z, n_z);
for i = 1:n_z
    for j = 1:n_z
        if j == 1
            transition_matrix(i, j) = normcdf((z_values(j) + z_interval / 2 - rho * z_values(i)) / sigma);
        elseif j == n_z
            transition_matrix(i, j) = 1 - normcdf((z_values(j) - z_interval / 2 - rho * z_values(i)) / sigma);
        else
            transition_matrix(i, j) = normcdf((z_values(j) + z_interval / 2 - rho * z_values(i)) / sigma) - ...
                normcdf((z_values(j) - z_interval / 2 - rho * z_values(i)) / sigma);
        end
    end
end

end
