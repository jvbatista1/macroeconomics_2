clear all; clc; close all;

% Parâmetros do Modelo
s = 2;              % Elasticidade da utilidade intertemporal
rho = 0.05;         % Taxa de impaciência
r = 0.03;           % Taxa de juros
z1 = 0.1;           % Produtividade do trabalho no estado 1
z2 = 0.2;           % Produtividade do trabalho no estado 2
z = [z1, z2];       % Vetor de produtividades do trabalho
la1 = 1.5;          % Parâmetro da função de utilidade
la2 = 1;            % Parâmetro da função de utilidade
la = [la1, la2];     % Vetor de parâmetros da função de utilidade
frisch = 1/2;       % Elasticidade da oferta de trabalho em relação ao salário real
w = 1;              % Salário real

% Configurações da grade de riqueza
I = 500;            % Número de pontos na grade
amin = -0.15;       % Riqueza mínima
amax = 3;           % Riqueza máxima
a = linspace(amin, amax, I)'; % Vetor de riqueza
da = (amax - amin)/(I - 1);   % Incremento na riqueza

% Inicialização de variáveis
aa = [a, a];           % Vetor de riqueza replicado para os dois estados
zz = ones(I, 1) * z;   % Vetor de produtividades replicado para cada ponto na grade

maxit = 20;            % Número máximo de iterações
crit = 1e-6;           % Critério de convergência
Delta = 1000;          % Parâmetro Delta para suavização
x0 = (w * z1)^(frisch*(1-s)/(1+s*frisch)); % Chute inicial para oferta de trabalho

% Inicialização de matrizes
dVf = zeros(I, 2);
dVb = zeros(I, 2);
c = zeros(I, 2);

% Matriz de transição para o estado exógeno
Aswitch = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];

% Configuração das opções para a função fzero
options = optimset('Display', 'off');

% Resolução do modelo   
tic;
for i=1:I
    % Encontrar oferta de trabalho para o estado 1
    params = [a(i),z1,w,r,s,frisch];
    myfun = @(l) lab_solve(l,params);
    [l01,fval,exitflag] = fzero(myfun,x0,options);

    % Encontrar oferta de trabalho para o estado 2
    params = [a(i),z2,w,r,s,frisch];
    myfun = @(l) lab_solve(l,params);
    [l02,fval,exitflag] = fzero(myfun,x0,options);

    % Armazenar ofertas de trabalho
    l0(i,:)=[l01,l02];
end
toc

% Calcular funções de valor iniciais
v0(:,1) = (w*z(1).*l0(1,1) + r.*a).^(1-s)/(1-s)/rho;
v0(:,2) = (w*z(2).*l0(1,2) + r.*a).^(1-s)/(1-s)/rho;
plot(a,v0(:,2))

% Iterações para resolver o modelo
lmin = l0(1,:);
lmax = l0(I,:);

v = v0;

%maxit = 1;
for n=1:maxit
    V = v;
    V_n(:,:,n)=V;
    % forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w*z.*lmax + r.*amax).^(-s); %state constraint boundary condition
    % backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w*z.*lmin + r.*amin).^(-s); %state constraint boundary condition
       
    %consumption and savings with forward difference
    cf = dVf.^(-1/s);
    lf = (dVf.*w.*zz).^frisch;
    ssf = w*zz.*lf + r.*aa - cf;
    %consumption and savings with backward difference
    cb = dVb.^(-1/s);
    lb = (dVb.*w.*zz).^frisch;
    ssb = w*zz.*lb + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = w*zz.*l0 + r.*aa;
    dV0 = c0.^(-s);
    
    Ib = ssb < 0; %negative drift --> backward difference
    If = (ssf > 0).*(1-Ib); %positive drift --> forward difference
    I0 = (1-If-Ib); %at steady state
    
    c = cf.*If + cb.*Ib + c0.*I0;
    l = lf.*If + lb.*Ib + l0.*I0;
    u = c.^(1-s)/(1-s) - l.^(1+1/frisch)/(1+1/frisch);
     
    %CONSTRUCT MATRIX
    X = -Ib.*ssb/da;
    Y = -If.*ssf/da + Ib.*ssb/da;
    Z = If.*ssf/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);

    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
    end
    
    B = (1/Delta + rho)*speye(2*I) - A;
    
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
      
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
toc;


adot = w.*zz.*l + r.*aa - c;

% Plotagem dos resultados
plot(a,adot,a,zeros(1,I),'--','LineWidth',2)
set(gca,'FontSize',14)
legend('s_1(a)','s_2(a)','Location','SouthWest')
grid
xlabel('Wealth, a')
ylabel('Saving, s_i(a)')
print -depsc labor_s.eps

plot(a,v,'LineWidth',2)
set(gca,'FontSize',14)
legend('v_1(a)','v_2(a)','Location','SouthWest')
grid
xlabel('Wealth, a')
ylabel('Value Function, v_i(a)')
print -depsc labor_v.eps

plot(a,l,'LineWidth',2)
set(gca,'FontSize',14)
legend('l_1(a)','l_2(a)','Location','SouthWest')
grid
xlabel('Wealth, a')
ylabel('Labor Supply, l_i(a)')
print -depsc labor_l.eps

plot(a,c,'LineWidth',2)
set(gca,'FontSize',14)
legend('c_1(a)','c_2(a)','Location','SouthWest')
grid
xlabel('Wealth, a')
ylabel('Consumption, c_i(a)')
print -depsc labor_c.eps