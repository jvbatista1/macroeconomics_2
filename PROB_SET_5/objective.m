function output = objective(x, mu_function, w_values, beta)
    % This function computes the objective function to be minimized
    % for the given parameters x = [lambda; mu]
    
    % Economic names for the parameters
    lambda = x(1);
    mu = x(2);
    
    % Calculate the stationary unemployment rate
    U = lambda / (lambda + mu);
    
    % Calculate the job destruction rate
    E2U = lambda;
    
    % Calculate the expected value of the wage offer
    expected_wage = mu * trapz(w_values .* lognpdf(w_values, 0, 1));
    
    % Compute the objective function
    output = mu_function(U) * trapz(max(expected_wage / (1 - beta * (1 - U)) + beta * U / (1 - beta * (1 - U)) * mu_function(U), beta * mu_function(U))) + (1 - mu_function(U)) * (beta * mu_function(U));
end