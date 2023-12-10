function equilibrium = find_equilibrium()
    % Parameters
    beta = 0.9;
    mu0 = 0.1; % You can adjust this value
    
    % Discretize the distribution of wage offers
    w_values = linspace(0, 20, 41);
    
    % Define the function for mu(U)
    mu_function = @(U) mu0 * (1 - sqrt(U));
    
    % Solve for equilibrium
    guess = [0.025; 0.50]; % Initial guess for lambda and mu
    
    % Set options for the minimization algorithm
    options = optimset('display', 'iter');
    
    % Call fminsearch to find the equilibrium
    equilibrium_params = fminsearch(@(x) objective(x, mu_function, w_values, beta), guess, options);
    
    % Display the equilibrium parameters
    lambda_eq = equilibrium_params(1);
    mu_eq = equilibrium_params(2);
    
    fprintf('Equilibrium Lambda: %.4f\n', lambda_eq);
    fprintf('Equilibrium Mu: %.4f\n', mu_eq);
    
    % Return the equilibrium parameters
    equilibrium.lambda = lambda_eq;
    equilibrium.mu = mu_eq;
end

