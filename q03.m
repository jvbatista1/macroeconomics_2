%%%%%%%%%%%%
% Macroeconomia II
% Lista I
% Questão 3
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%

%% (a)
% Declarando a função f(z) = e^z
f = @(z) (z >= 0 & z <= 1) .* exp(z);

% Definindo a ordem e os pontos de Chebyshev
n = 2;
m = 3;

% Calculando os nós de Chebyshev manualmente
cheb_nodes = cos((2*(1:m)-1)*pi/(2*m));

% Avaliando a função nos pontos de Chebyshev
cheb_values = f(cheb_nodes);

% Ajustando um polinômio usando polyfit
cheb_coeffs = polyfit(cheb_nodes, cheb_values, n);

% Criando uma função a partir dos coeficientes do polinômio
cheb_approx = @(z) polyval(cheb_coeffs, z);

% Criando um vetor de pontos para o gráfico
z_values = linspace(0, 1, 1000);

% Calculando os valores da função original e da aproximação nos pontos
f_values = f(z_values);
approx_values = cheb_approx(z_values);

% Plotando a função original e a aproximação
figure;
plot(z_values, f_values, 'LineWidth', 2, 'DisplayName', 'f(z) = e^z');
hold on;
plot(z_values, approx_values, '--', 'LineWidth', 2, 'DisplayName', 'Aproximação de Chebyshev');
scatter(cheb_nodes, cheb_values, 100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Pontos de Chebyshev');
legend('show');
xlabel('z');
ylabel('f(z)');
title('Interpolação de Chebyshev');
grid on;
hold off;

%% (b)
% Definindo a função original
f = @(z) exp(z);

% Definindo a ordem e os pontos de Chebyshev
n = 2;
m1 = 3; % pontos de Chebyshev para a primeira interpolação
m2 = 10; % pontos de Chebyshev para a segunda interpolação

% Calculando os nós de Chebyshev manualmente
cheb_nodes1 = cos((2*(1:m1)-1)*pi/(2*m1));
cheb_nodes2 = cos((2*(1:m2)-1)*pi/(2*m2));

% Avaliando a função nos pontos de Chebyshev
cheb_values1 = f(cheb_nodes1);
cheb_values2 = f(cheb_nodes2);

% Ajustando polinômios usando polyfit
cheb_coeffs1 = polyfit(cheb_nodes1, cheb_values1, n);
cheb_coeffs2 = polyfit(cheb_nodes2, cheb_values2, n);

% Criando funções a partir dos coeficientes dos polinômios
cheb_approx1 = @(z) polyval(cheb_coeffs1, z);
cheb_approx2 = @(z) polyval(cheb_coeffs2, z);

% Criando um vetor de pontos para o gráfico
z_values = linspace(0, 1, 1000);

% Calculando os valores das funções originais e das aproximações nos pontos
f_values = f(z_values);
approx_values1 = cheb_approx1(z_values);
approx_values2 = cheb_approx2(z_values);

% Plotando as funções originais e as aproximações
figure;
plot(z_values, f_values, 'LineWidth', 2, 'DisplayName', 'f(z) = e^z');
hold on;
plot(z_values, approx_values1, '--', 'LineWidth', 2, 'DisplayName', 'Aproximação 1 de Chebyshev');
plot(z_values, approx_values2, '-.', 'LineWidth', 2, 'DisplayName', 'Aproximação 2 de Chebyshev');
scatter(cheb_nodes1, cheb_values1, 100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Pontos 1 de Chebyshev');
scatter(cheb_nodes2, cheb_values2, 100, 'b', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Pontos 2 de Chebyshev');
legend('show');
xlabel('z');
ylabel('f(z)');
title('Interpolação de Chebyshev');
grid on;
hold off;

%% (c)
