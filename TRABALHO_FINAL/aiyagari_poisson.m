%Optimized for speed by SeHyoun Ahn

clear all; clc; close all;

tic;

% Parâmetros do Modelo
ga = 2;% Coeficiente de aversão ao risco
rho = 0.05;% Taxa de desconto
d = 0.05;% Taxa de depreciação do capital
al = 1/3;       % Parâmetro de elasticidade da função de produção
Aprod = 0.1;    % Produtividade total dos fatores
z1 = 1;         % Nível de produtividade para o tipo 1 de agente
z2 = 2 * z1;    % Nível de produtividade para o tipo 2 de agente
z = [z1, z2];   % Vetor de níveis de produtividade
la1 = 1/3;      % Fração da população para o tipo 1 de agente
la2 = 1/3;      % Fração da população para o tipo 2 de agente
la = [la1, la2]; % Vetor de frações populacionais
z_ave = (z1 * la2 + z2 * la1) / (la1 + la2); % Média ponderada da produtividade

% Configuração da Grade de Ativos
I = 1000;       % Número de pontos na grade de ativos
amin = 0;       % Limite inferior da grade de ativos
amax = 20;      % Limite superior da grade de ativos

% Configuração Inicial da Grade e Matrizes
a = linspace(amin, amax, I)'; % Vetor de ativos
da = (amax - amin) / (I - 1); % Tamanho do passo da grade
aa = [a, a];                  % Matriz de ativos repetida para cada tipo de agente
zz = ones(I, 1) * z;          % Matriz de produtividade para cada ponto na grade


% Parâmetros do Loop de Resolução do Modelo
maxit = 100;   % Número máximo de iterações
crit = 1e-6;   % Critério de convergência para a função de valor
Delta = 1000;  % Parâmetro de iteração

dVf = zeros(I, 2); % Inicialização da matriz forward difference
dVb = zeros(I, 2); % Inicialização da matriz backward difference
c = zeros(I, 2);   % Inicialização da matriz de consumo

Aswitch = [-speye(I) * la(1), speye(I) * la(1); speye(I) * la(2), -speye(I) * la(2)];

Ir = 40;       % Número de iterações no loop externo
crit_S = 1e-5;  % Critério de convergência externo

% Configuração Inicial da Taxa de Juros
rmax = 0.049;
r = 0.04;
w = 0.05;

r0 = 0.03;
rmin = 0.01;   % Limite inferior para a taxa de juros
rmax = 0.99 * rho; % Limite superior para a taxa de juros

% Loop para a Resolução do Modelo
for ir=1:Ir
    r_r(ir)=r;
    rmin_r(ir)=rmin;
    rmax_r(ir)=rmax;

    % Cálculo da Demanda de Capital e Salário
    KD(ir) = (al*Aprod/(r + d))^(1/(1-al))*z_ave;
    w = (1-al)*Aprod*KD(ir).^al*z_ave^(-al);

    % Verificação da Restrição de Empréstimo
    if w*z(1) + r*amin < 0
        disp('CAREFUL: borrowing constraint too loose')
    end

    % Inicialização da Função de Valor
    v0(:,1) = (w*z(1) + r.*a).^(1-ga)/(1-ga)/rho;
    v0(:,2) = (w*z(2) + r.*a).^(1-ga)/(1-ga)/rho;

    % Utiliza a iteração anterior se ir > 1
    if ir > 1
        v0 = V_r(:,:,ir-1);
    end

    v = v0;

    % Iteração de Valor
    for n=1:maxit
        V = v;
        V_n(:,:,n)=V;
        % Diferença para Frente
        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        dVf(I,:) = (w*z + r.*amax).^(-ga); % condição de estado a<=amax
        % Diferença para Trás
        dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
        dVb(1,:) = (w*z + r.*amin).^(-ga); % condição de estado a>=amin

        % Consumo e Poupança com Diferença para Frente
        cf = dVf.^(-1/ga);
        ssf = w*zz + r.*aa - cf;
        % Consumo e Poupança com Diferença para Trás
        cb = dVb.^(-1/ga);
        ssb = w*zz + r.*aa - cb;
        % Consumo e Derivada da Função de Valor em Estado Estacionário
        c0 = w*zz + r.*aa;

        % Escolha de Diferença para Frente ou para Trás com base na direção do drift    
        If = ssf > 0; % drift positivo --> diferença para frente
        Ib = ssb < 0; % drift negativo --> diferença para trás
        I0 = (1-If-Ib); % em estado estacionário

        c = cf.*If + cb.*Ib + c0.*I0;
        u = c.^(1-ga)/(1-ga);

        % Construção da Matriz Tridiagonal
        X = -min(ssb,0)/da;
        Y = -max(ssf,0)/da + min(ssb,0)/da;
        Z = max(ssf,0)/da;

        A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;

        % Verifica se a matriz de transição é adequada
        if max(abs(sum(A,2)))>10^(-9)
           disp('Improper Transition Matrix')
           %break
        end

        B = (1/Delta + rho)*speye(2*I) - A;

        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];

        b = u_stacked + V_stacked/Delta;
        V_stacked = B\b; % Solução do Sistema de Equações

        V = [V_stacked(1:I),V_stacked(I+1:2*I)];

        Vchange = V - v;
        v = V;

        dist(n) = max(max(abs(Vchange)));
        if dist(n)<crit
            disp('Função de Valor Convergida, Iteração = ')
            disp(n)
            break
        end
    end

    % FOKKER-PLANCK EQUATION
    AT = A';
    b = zeros(2*I,1);

    % Corrige um valor para evitar singularidade na matriz
    i_fix = 1;
    b(i_fix) = 0.1;
    row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
    AT(i_fix,:) = row;

    % Resolve o sistema linear
    gg = AT\b;
    g_sum = gg'*ones(2*I,1)*da;
    gg = gg./g_sum;

    g = [gg(1:I),gg(I+1:2*I)];

    check1 = g(:,1)'*ones(I,1)*da;
    check2 = g(:,2)'*ones(I,1)*da;

    g_r(:,:,ir) = g;
    adot(:,:,ir) = w*zz + r.*aa - c;
    V_r(:,:,ir) = V;

    KS(ir) = g(:,1)'*a*da + g(:,2)'*a*da;
    S(ir) = KS(ir) - KD(ir);

    % Atualização da Taxa de Juros
    if S(ir)>crit_S
        disp('Excesso de Oferta')
        rmax = r;
        r = 0.5*(r+rmin);
    elseif S(ir)<-crit_S
        disp('Excesso de Demanda')
        rmin = r;
        r = 0.5*(r+rmax);
    elseif abs(S(ir))<crit_S
        display('Equilíbrio Encontrado, Taxa de Juros =')
        disp(r)
        break
    end
end

% Plotagem dos Resultados
amax1 = 5;
amin1 = amin - 0.1;

% Figura 1: Funções de Poupança
figure(1)
h1 = plot(a, adot(:,1,ir), 'b', a, adot(:,2,ir), 'r', linspace(amin1,amax1,I), zeros(1,I), 'k--', 'LineWidth', 2);
legend(h1, 's_1(a)', 's_2(a)', 'Location', 'NorthEast');
text(-0.155, -0.105, '$\underline{a}$', 'FontSize', 16, 'interpreter', 'latex');
line([amin amin], [-0.1 0.08], 'Color', 'Black', 'LineStyle', '--');
xlabel('Riqueza, $a$', 'interpreter', 'latex');
ylabel('Poupança, $s_i(a)$', 'interpreter', 'latex');
xlim([amin1 amax1]);
ylim([-0.03 0.05]);
set(gca, 'FontSize', 16);

% Figura 2: Funções de Densidade
figure(2)
h1 = plot(a, g_r(:,1,ir), 'b', a, g_r(:,2,ir), 'r', 'LineWidth', 2);
legend(h1, 'g_1(a)', 'g_2(a)');
text(-0.155, -0.12, '$\underline{a}$', 'FontSize', 16, 'interpreter', 'latex');
line([amin amin], [0 max(max(g_r(:,:,ir)))], 'Color', 'Black', 'LineStyle', '--');
xlabel('Riqueza, $a$', 'interpreter', 'latex');
ylabel('Densidades, $g_i(a)$', 'interpreter', 'latex');
xlim([amin1 amax1]);
%ylim([0 0.5])
set(gca, 'FontSize', 16);

toc;
