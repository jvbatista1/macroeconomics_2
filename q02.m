%%%%%%%%%%%%
% Macroeconomia II
% Lista I
% Questão 2
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%

%% (a.1)
% Defining c
c = 1;

% Defining d(x, c)
d = @(x) x.^(-5) - x.^(-3) - c;

% Different values of x
x_values = linspace(0.6, 10, 100);  
result = d(x_values);

% Plot
figure;
plot(x_values, result, '-o', 'LineWidth', 2);
xlabel('x');
ylabel('d(x)');
title('Gráfico de d(x) = x^{-5} - x^{-3} - c');
grid on;

% indication of c value
legend(['c = ' num2str(c)], 'Location', 'Best');

%% (a.2)
% Definig the function and its first derivative
d = @(x, c) x.^(-5) - x.^(-3) - c;
d_prime = @(x) 5*x.^(-6) - 3*x.^(-4);

% X initial value
x0 = 0.5;

% Precision criteria
tolerance = 1e-8;

% Max iterations
max_iterations = 1000;

% Newton-Raphson iteration
x = x0;
iteration = 0;
while abs(d(x, c)) > tolerance && iteration < max_iterations
    x = x - d(x, c) / d_prime(x);
    iteration = iteration + 1;
end

root_approximation = x;

disp(['Root approximation: ' num2str(root_approximation)]);

%% (b)
d = @(x, c) x.^(-5) - x.^(-3) - c;
d_prime = @(x) 5*x.^(-6) - 3*x.^(-4);

% Definig c values grid
c_values = linspace(1, 8, 10);

x0 = 0.5;
tolerance = 1e-8;
max_iterations = 1000;

% Vector to store different root approx
root_approximations = zeros(size(c_values));

% Iteration to different c values
for i = 1:length(c_values)
    % Initialize for each iteration
    x = x0;
    iteration = 0;

    % NR function
    while abs(d(x, c_values(i))) > tolerance && iteration < max_iterations
        x = x - d(x, c_values(i)) / d_prime(x);
        iteration = iteration + 1;
    end

    % Storing
    root_approximations(i) = x;
    
    % Displaying
    disp(['For c = ' num2str(c_values(i)) ', root: ' num2str(root_approximations(i))]);
end


%% (c)
d = @(x, c) x.^(-5) - x.^(-3) - c;
d_prime = @(x) 5*x.^(-6) - 3*x.^(-4);

% New c values
c_values = linspace(1, 10, 1000);
x0 = 0.5;
tolerance = 1e-8;
max_iterations = 1000;
root_approximations = zeros(size(c_values));

for i = 1:length(c_values)
    x = x0;
    iteration = 0;
    while abs(d(x, c_values(i))) > tolerance && iteration < max_iterations
        x = x - d(x, c_values(i)) / d_prime(x);
        iteration = iteration + 1;
    end
    root_approximations(i) = x;
end

% Cubic spline
c_spline = linspace(1, 10, 1000);
x_spline = spline(c_values, root_approximations, c_spline);

figure;
plot(c_values, root_approximations, '-o', 'LineWidth', 2, 'DisplayName', 'Aproximação de x*');
hold on;
plot(c_spline, x_spline, '--', 'LineWidth', 2, 'DisplayName', 'Spline Cúbica');
xlabel('c');
ylabel('Aproximação de x*');
title('Aproximação de x* para diferentes valores de c usando o Método de Newton-Raphson com Spline Cúbica');
grid on;
legend('show');