% Limpeza do ambiente
clear all; clc; close all;

% Parâmetros do Modelo
ga = 2;              % Coeficiente de aversão ao risco
rho = 0.05;          % Taxa de impaciência intertemporal
d = 0.05;            % Taxa de depreciação do capital
al = 1/3;            % Participação do capital na produção
Aprod = 0.1;         % Nível inicial de produtividade total dos fatores
z1 = 1;              % Estado de produtividade 1
z2 = 2*z1;           % Estado de produtividade 2
z = [z1, z2];        % Vetor de estados de produtividade
la1 = 1/3;           % Probabilidade do estado de produtividade 1
la2 = 1/3;           % Probabilidade do estado de produtividade 2
la = [la1, la2];      % Vetor de probabilidades
z_ave = (z1*la2 + z2*la1)/(la1 + la2);  % Média ponderada dos estados de produtividade
% Adicione as taxas de imposto heterogêneas para cada estado de produtividade
tax_rates = [0.1, 0.2];
tau_z1 = tax_rates(1);
tau_z2 = tax_rates(2);

% Parâmetros Temporais e de Discretização
T = 200;             % Tempo total
N = 400;             % Número de intervalos de tempo
dt = T/N;            % Passo de tempo
time = (0:N-1)*dt;    % Vetor de tempo
max_price_it = 300;   % Número máximo de iterações para encontrar preços de equilíbrio
convergence_criterion = 10^(-5);  % Critério de convergência para a transição dinâmica
relax = 0.1;          % Parâmetro de relaxamento

% Construção da Sequência de Produtividade (Processo de Poisson)
corr = 0.8;          % Correlação entre os estados
nu = 1 - corr;       % Parâmetro relacionado à construção do processo estocástico
Aprod_t = zeros(N,1); % Sequência de produtividade ao longo do tempo
Aprod_t(1) = 0.97*Aprod;
for n = 1:N-1
    Aprod_t(n+1) = dt*nu*(Aprod - Aprod_t(n)) + Aprod_t(n);
end

% Plot da Sequência de Produtividade ao Longo do Tempo
plot(time, Aprod_t)
xlim([0 40])
xlabel('Tempo')
ylabel('Produtividade')
title('Sequência de Produtividade ao Longo do Tempo')

% Parâmetros para a Grade de Capital
I = 1000;            % Número de pontos na grade de acumulação de capital
amin = -0.8;         % Valor mínimo para o capital
amax = 20;           % Valor máximo para o capital
a = linspace(amin, amax, I)'; % Vetor de capital
da = (amax - amin)/(I-1);      % Incremento na grade de capital
aa = [a, a];        % Vetor de capital repetido para cada estado de produtividade
zz = ones(I, 1)*z;   % Vetor de estados de produtividade repetido para cada ponto na grade

% Inicialização de Variáveis para a Iteração da Função Valor
maxit = 100;           % Número máximo de iterações
crit = 10^(-5);        % Critério de convergência para a função valor
Delta = 1000;          % Parâmetro para iteração da função valor
dVf = zeros(I, 2);     % Forward differences da função valor
dVb = zeros(I, 2);     % Backward differences da função valor
c = zeros(I, 2);       % Consumo
Aswitch = [-speye(I)*la(1), speye(I)*la(1); speye(I)*la(2), -speye(I)*la(2)]; % Matriz de transição entre estados de produtividade

% Parâmetros para Encontrar o Equilíbrio Estacionário
Ir = 40;               % Número de iterações para encontrar o equilíbrio estacionário
crit_S = 10^(-5);      % Critério de convergência para a oferta e demanda de capital
rmax = 0.049;          % Taxa de juros máxima
r = 0.04;              % Taxa de juros inicial
w = 0.05;              % Salário inicial

% Inicialização das Condições Iniciais da Função Valor
v0(:, 1) = (w*z(1) + r*a).^(1-ga)/(1-ga)/rho;
v0(:, 2) = (w*z(2) + r*a).^(1-ga)/(1-ga)/rho;

r0 = 0.03;             % Taxa de juros inicial (usado para verificação de convergência)
rmin = 0.01;           % Taxa de juros mínima
rmax = 0.99*rho;       % Taxa de juros máxima (usado para verificação de convergência)


%%%%%%%%%%%%%%%%
% STEADY STATE %
%%%%%%%%%%%%%%%%

% Iteração para encontrar o equilíbrio estacionário
for ir = 1:Ir
    r_r(ir) = r;
    rmin_r(ir) = rmin;
    rmax_r(ir) = rmax;

    % Cálculo do capital demandado e salário
    KD(ir) = (al*Aprod/(r + d))^(1/(1-al))*z_ave;
    w = (1-al)*Aprod*KD(ir).^al*z_ave^(-al);

    % Inicialização da função valor se ir > 1
    if ir > 1
        v0 = V_r(:, :, ir-1);
    end

    v = v0;

    % Iteração da função valor
    for n = 1:maxit
        V = v;
        V_n(:, :, n) = V;

        % Forward difference
        dVf(1:I-1, :) = (V(2:I, :) - V(1:I-1, :))/da;
        dVf(I, :) = (w*z + r*amax).^(-ga);

        % Backward difference
        dVb(2:I, :) = (V(2:I, :) - V(1:I-1, :))/da;
        dVb(1, :) = (w*z + r*amin).^(-ga);



        % % Consumo e poupança com forward difference
        % cf = dVf.^(-1/ga);
        % ssf = w*zz + r*aa - cf;
        % 
        % % Consumo e poupança com backward difference
        % cb = dVb.^(-1/ga);
        % ssb = w*zz + r*aa - cb;
        %
        % % Consumo e derivada da função valor no estado estacionário
        % c0 = w*zz + r*aa;
        % dV0 = c0.^(-ga);
        % 
        % % Escolha do método de diferenças baseado no sinal do drift
        % If = ssf > 0; % Forward difference
        % Ib = ssb < 0; % Backward difference
        % I0 = (1 - If - Ib); % Estado estacionário
        % 
        % % Diferenças baseadas no método upwind
        % dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0;
        % c = dV_Upwind.^(-1/ga);
        % u = c.^(1-ga)/(1-ga);

        % Consumo e poupança com forward difference (estado de produtividade z1)
        cf_z1 = (1 - tau_z1)*dVf(:,1).^(-1/ga);
        cb_z1 = (1 - tau_z1)*dVb(:,1).^(-1/ga);
    
        ssf_z1 = w*z(1) + r*aa - cf_z1;
        ssb_z1 = w*z(1) + r*aa - cb_z1;

        % Consumo e poupança com backward difference (estado de produtividade z1)
        cf_z2 = (1 - tau_z2)*dVf(:,2).^(-1/ga);
        cb_z2 = (1 - tau_z2)*dVb(:,2).^(-1/ga);
    
        ssf_z2 = w*z(2) + r*aa - cf_z2;
        ssb_z2 = w*z(2) + r*aa - cb_z2;

        % Consumo e derivada da função valor no estado estacionário (estado de produtividade z1)
        c0_z1 = w*z(1) + r*aa;
        dV0_z1 = c0_z1.^(-ga);
        
        % Consumo e derivada da função valor no estado estacionário (estado de produtividade z2)
        c0_z2 = w*z(2) + r*aa;
        dV0_z2 = c0_z2.^(-ga);
    
        % Escolha do método de diferenças baseado no sinal do drift (estado de produtividade z1)
        If_z1 = ssf_z1 > tau_z1; % Forward difference
        Ib_z1 = ssb_z1 < tau_z1; % Backward difference
        I0_z1 = (1 - If_z1 - Ib_z1); % Estado estacionário
    
        % Diferenças baseadas no método upwind (estado de produtividade z1)
        dV_Upwind_z1 = dVf.*If_z1 + dVb.*Ib_z1 + dV0_z1.*I0_z1;
        c_z1 = (dV_Upwind_z1 - tau_z1*w*z(1)).^(-1/ga);
        u_z1 = c_z1.^(1-ga)/(1-ga);
    
        % Escolha do método de diferenças baseado no sinal do drift (estado de produtividade z2)
        If_z2 = ssf_z2 > tau_z2; % Forward difference
        Ib_z2 = ssb_z2 < tau_z2; % Backward difference
        I0_z2 = (1 - If_z2 - Ib_z2); % Estado estacionário
    
        % Diferenças baseadas no método upwind (estado de produtividade z2)
        dV_Upwind_z2 = dVf.*If_z2 + dVb.*Ib_z2 + dV0_z2.*I0_z2;
        c_z2 = (dV_Upwind_z2 - tau_z2*w*z(2)).^(-1/ga);
        u_z2 = c_z2.^(1-ga)/(1-ga);

        % Construção da matriz para z1
        X_z1 = -min(ssb_z1, 0)/da;
        Y_z1 = -max(ssf_z1, 0)/da + min(ssb_z1, 0)/da;
        Z_z1 = max(ssf_z1, 0)/da;
        
        %A1_z1 = spdiags(Y_z1, 0, I, I) + spdiags(X_z1(2:I), -1, I, I) + spdiags([0; Z_z1(1:I-1)], 1, I, I);
        A1_z1 = spdiags([Y_z1, X_z1, Z_z1], [0, -1, 1], I, I);
        
        
        % Construção da matriz para z2
        X_z2 = -min(ssb_z2, 0)/da;
        Y_z2 = -max(ssf_z2, 0)/da + min(ssb_z2, 0)/da;
        Z_z2 = max(ssf_z2, 0)/da;
        
        %A2_z2 = spdiags(Y_z2, 0, I, I) + spdiags(X_z2(2:I), -1, I, I) + spdiags([0; Z_z2(1:I-1)], 1, I, I);
        A2_z2 = spdiags([Y_z2, X_z2, Z_z2], [0, -1, 1], I, I);
        
        % Montagem da matriz A
        A = [A1_z1, sparse(I, I); sparse(I, I), A2_z2] + Aswitch;
        
        % Verificação da matriz de transição
        if max(abs(sum(A, 2))) > 10^(-9)
            disp('Improper Transition Matrix')
            break
        end
        
        % Construção da matriz do sistema
        B = (1/Delta + rho)*speye(2*I) - A;
        
        u_stacked = [u_z1; u_z2];
        V_stacked = [V(:, 1); V(:, 2)];
        
        b = u_stacked + V_stacked/Delta;
        V_stacked = B\b; % Solução do sistema de equações
        
        V = [V_stacked(1:I), V_stacked(I+1:2*I)];
        
        % Verificação da convergência
        Vchange = V - v;
        v = V;
        
        dist(n) = max(max(abs(Vchange)));
        
        if dist(n) < crit
            disp('Função Valor Convergida, Iteração = ')
            disp(n)
            break
        end
    end
    toc;

    % Construção da matriz transposta de A para a equação de Fokker-Planck
    AT = A';
    b = zeros(2*I, 1);

    % Necessário fixar um valor para evitar matriz singular
    i_fix = 1;
    b(i_fix) = 0.1;
    row = [zeros(1, i_fix-1), 1, zeros(1, 2*I-i_fix)];
    AT(i_fix, :) = row;

    % Resolução do sistema linear
    gg = AT\b;
    g_sum = gg'*ones(2*I, 1)*da;
    gg = gg./g_sum;

    g = [gg(1:I), gg(I+1:2*I)];

    % Verificação de algumas propriedades da distribuição
    check1 = g(:, 1)'*ones(I, 1)*da;
    check2 = g(:, 2)'*ones(I, 1)*da;

    g_r(:, :, ir) = g;
    adot(:, :, ir) = w*zz + r*aa - c;
    V_r(:, :, ir) = V;

    KS(ir) = g(:,1)'*a*da + g(:,2)'*a*da;
    S(ir) = KS(ir) - KD(ir);

    % Atualização da taxa de juros
    if S(ir) > crit_S
        disp('Excesso de Oferta')
        rmax = r;
        r = 0.5*(r + rmin);
    elseif S(ir) < -crit_S
        disp('Excesso de Demanda')
        rmin = r;
        r = 0.5*(r + rmax);
    elseif abs(S(ir)) < crit_S
        disp('Equilíbrio Encontrado, Taxa de Juros =')
        disp(r)
        break
    end

end