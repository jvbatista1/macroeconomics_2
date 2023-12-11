function r_equilibrium = find_equilibrium_interest_rate(beta, gamma, w, phi, r_low, r_high)
    epsilon = 1e-6; % Critério de convergência

    while abs(r_high - r_low) > epsilon
        r_guess = (r_low + r_high) / 2;

        % Calcular a demanda de ativos
        asset_demand = calculate_asset_demand(beta, gamma, w, phi, r_guess);

        % Atualizar r_high ou r_low com base na demanda
        if asset_demand > 0
            r_high = r_guess;
        else
            r_low = r_guess;
        end
    end

    % O equilíbrio está no meio entre r_low e r_high
    r_equilibrium = (r_low + r_high) / 2;
end
